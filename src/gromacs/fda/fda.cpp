/*
 * fda.cpp
 *
 *  Created on: Oct 31, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <sstream>
#include "fda.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "pf_array.h"
#include "pf_interactions.h"
#include "Utilities.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

static const real HALF    = 1.0 / 2.0;
static const real THIRD   = 1.0 / 3.0;
static const real QUARTER = 0.25;

FDA::FDA(FDASettings const& fda_settings)
 : fda_settings(fda_settings),
   atom_based(fda_settings.atom_based_result_type,
		      fda_settings.one_pair,
			  fda_settings.syslen_atoms,
		      fda_settings.atom_based_result_filename,
			  fda_settings.no_end_zeros,
			  fda_settings.v2s),
   residue_based(fda_settings.residue_based_result_type,
		         fda_settings.one_pair,
				 fda_settings.syslen_residues,
				 fda_settings.residue_based_result_filename,
		         fda_settings.no_end_zeros,
				 fda_settings.v2s),
   time_averaging_steps(0),
   time_averaging_com(nullptr),
   nsteps(0)
{
  if (atom_based.VS_mode()) atom_vir.resize(fda_settings.syslen_atoms);
}

void FDA::add_bonded_nocheck(int i, int j, int type, rvec force)
{
  // the calling functions will not have i == j, but there is not such guarantee for ri and rj;
  // and it makes no sense to look at the interaction of a residue to itself
  if (residue_based.PF_or_PS_mode()) {
    int ri = fda_settings.get_atom2residue(i);
    int rj = fda_settings.get_atom2residue(j);
    rvec force_residue;
    if (ri != rj) {
      switch(fda_settings.one_pair) {
        case OnePair::DETAILED:
          if (ri > rj) {
            int_swap(&ri, &rj);
            clear_rvec(force_residue);
            rvec_dec(force_residue, force);
            residue_based.distributed_forces.detailed[ri][rj][type] += force_residue;
          } else {
            residue_based.distributed_forces.detailed[ri][rj][type] += force;
          }
          break;
        case OnePair::SUMMED:
          if (ri > rj) {
            int_swap(&ri, &rj);
            clear_rvec(force_residue);
            rvec_dec(force_residue, force);
            residue_based.distributed_forces.summed[ri][rj] += force_residue;
          } else {
            residue_based.distributed_forces.summed[ri][rj] += force;
          }
          break;
        case OnePair::INVALID:
          break;
      }
    }
  }

  if (atom_based.PF_or_PS_mode()) {
    if (i > j) {
      int_swap(&i, &j);
      rvec_opp(force);
    }
    switch(fda_settings.one_pair) {
      case OnePair::DETAILED:
        atom_based.distributed_forces.detailed[i][j][type] += force;
        break;
      case OnePair::SUMMED:
        atom_based.distributed_forces.summed[i][j] += force;
        break;
      case OnePair::INVALID:
        break;
    }
  }
}

void FDA::add_bonded(int i, int j, int type, rvec force)
{
  // leave early if the interaction is not interesting
  if (!(fda_settings.type & type)) return;
  if (!fda_settings.atoms_in_groups(i, j)) return;

  add_bonded_nocheck(i, j, type, force);
}

void FDA::add_nonbonded_single(int i, int j, int type, real force, real dx, real dy, real dz)
{
  // leave early if the interaction is not interesting
  if (!(fda_settings.type & type)) return;
  if (!fda_settings.atoms_in_groups(i, j)) return;

  rvec force_v;
  force_v[0] = force * dx;
  force_v[1] = force * dy;
  force_v[2] = force * dz;
  add_bonded_nocheck(i, j, type, force_v);
}

void FDA::add_angle(int ai, int aj, int ak, rvec f_i, rvec f_j, rvec f_k)
{
	rvec uf_i, uf_j, uf_k, f_j_i, f_j_k, f_i_k;
	//rvec wf, urik, rik;
	real nf_j_i, nf_j_k;
	//fprintf(stderr, "i=%d, j=%d, k=%d\n", ai+1, aj+1, ak+1);
	//fprintf(stderr, "f_i=%8f %8f %8f\n", f_i[0], f_i[1], f_i[2]);
	//fprintf(stderr, "f_j=%8f %8f %8f\n", f_j[0], f_j[1], f_j[2]);
	//fprintf(stderr, "f_k=%8f %8f %8f\n", f_k[0], f_k[1], f_k[2]);
	//rvec_sub(f_i, f_k, wf);
	//fprintf(stderr, "wf=%8f %8f %8f\n", wf[0], wf[1], wf[2]);
	//rvec_sub(x[ai], x[ak], rik);
	//unitv(rik, urik);
	unitv(f_i, uf_i);
	unitv(f_j, uf_j);
	unitv(f_k, uf_k);
	nf_j_i = - norm(f_i) * iprod(uf_i, uf_j);
	nf_j_k = - norm(f_k) * iprod(uf_k, uf_j);
	svmul(nf_j_i, uf_j, f_j_i);
	svmul(nf_j_k, uf_j, f_j_k);
	rvec_add(f_j_i, f_i, f_i_k);
	/*rvec_add(f_j_k, f_k, f_k_i);*/
	//fprintf(stderr, "f_j_i=%8f %8f %8f f_j_k=%8f %8f %8f, norm(f_j_i)=%8f, norm(f_j_k)=%8f\n", f_j_i[0], f_j_i[1], f_j_i[2], f_j_k[0], f_j_k[1], f_j_k[2], norm(f_j_i), norm(f_j_k));
	//fprintf(stderr, "f_i_k=%8f %8f %8f f_k_i=%8f %8f %8f\n", f_i_k[0], f_i_k[1], f_i_k[2], f_k_i[0], f_k_i[1], f_k_i[2]);
	//fprintf(stderr, "angle(urik, f_i_k)=%8f\n", acos(cos_angle(urik, f_i_k))*RAD2DEG);
	add_bonded(aj, ai, PF_INTER_ANGLE, f_j_i);
	add_bonded(ai, ak, PF_INTER_ANGLE, f_i_k);
	add_bonded(aj, ak, PF_INTER_ANGLE, f_j_k);
}

void FDA::add_dihedral(int i, int j, int k, int l, rvec f_i, rvec f_j, rvec f_k, rvec f_l)
{
    rvec f_mj, f_mk;
    rvec f_ipl, f_jpk, f_jpk_i, f_jpk_l, f_j_ipl, f_jpk_ipl;
    rvec f_j_i, f_k_i, f_j_l, f_k_l, f_j_k, f_l_i;
    rvec f_jpk_c_f_j, f_jpk_c_f_k;
    rvec uf_jpk, uf_j, uf_k;
    real cos_a, sin_a, cos_b, sin_b, sinacosbpsinbcosa;
    real nf_ipl, nf_jpk, nf_j, nf_k, nf_j_i, nf_j_l, nf_k_i, nf_k_l, nf_jpkxnf_j, nf_jpkxnf_k, nf_jpk_i, nf_jpk_l;

    /* below computation can sometimes return before finishing to avoid division with very small numbers;
     * this situation can occur f.e. when all f_i, f_j, f_k and f_l are (almost) zero;
     * in this case there is no call to pf_atom_add_bonded, no pairwise forces are recorded (which is different from recording zero forces!)
     */
    //fprintf(stderr, "dihedral: i=%d, j=%d, k=%d, l=%d\n", i, j, k, l);
    /*print_vec("f_i", f_i);
    print_vec("f_j", f_j);
    print_vec("f_k", f_k);
    print_vec("f_l", f_l);*/
    /* below computation needs -f_j and -f_k */
    clear_rvec(f_mj);
    rvec_sub(f_mj, f_j, f_mj);
    clear_rvec(f_mk);
    rvec_sub(f_mk, f_k, f_mk);
    //fprintf(stderr, "fi=%8f %8f %8f, fmj=%8f %8f %8f, fmk=%8f %8f %8f, fl=%8f %8f %8f\n", f_i[0], f_i[1], f_i[2], f_mj[0], f_mj[1], f_mj[2], f_mk[0], f_mk[1], f_mk[2], f_l[0], f_l[1], f_l[2]);
    /*print_vec("f_i", f_i);
    print_vec("f_j", f_mj);
    print_vec("f_k", f_mk);
    print_vec("f_l", f_l);*/
    rvec_add(f_i, f_l, f_ipl);
    rvec_add(f_mj, f_mk, f_jpk);
    unitv(f_jpk, uf_jpk);
    unitv(f_mj, uf_j);
    unitv(f_mk, uf_k);
    /*print_vec("uf_j", uf_j);
    print_vec("uf_k", uf_k);*/
    /* project f_i and f_l on f_jpk = f_j + f_k */
    nf_ipl = norm(f_ipl);
    /*fprintf(stderr, "nf_ipl=%f\n", nf_ipl);*/
    if (nf_ipl < GMX_FLOAT_EPS)
      return;
    svmul(iprod(f_i, f_ipl) / nf_ipl, uf_jpk, f_jpk_i);
    svmul(iprod(f_l, f_ipl) / nf_ipl, uf_jpk, f_jpk_l);
    rvec_add(f_jpk_i, f_jpk_l, f_jpk_ipl);
    /*print_vec("f_jpk_i", f_jpk_i);
    print_vec("f_jpk_l", f_jpk_l);
    print_vec("f_jpk_ipl", f_jpk_ipl);*/
    //fprintf(stderr, "f_jpk=%8f %8f %8f, f_jpk_ipl=%8f %8f %8f\n", f_jpk[0], f_jpk[1], f_jpk[2], f_jpk_ipl[0], f_jpk_ipl[1], f_jpk_ipl[2]);
    /* decompose f_jpk_i in 2 forces in directions of f_j and f_k; decompose f_jpk_l in 2 forces in directions of f_j and f_k */
    nf_jpk = norm(f_jpk);
    nf_j = norm(f_mj);
    nf_k = norm(f_mk);
    /* a = angle between f_jpk and -f_j, b = angle between f_jpk and -f_k */
    /* obtain cos from dot product and sin from cross product */
    cprod(f_jpk, f_mj, f_jpk_c_f_j);
    cprod(f_jpk, f_mk, f_jpk_c_f_k);
    nf_jpkxnf_j = nf_jpk * nf_j;
    nf_jpkxnf_k = nf_jpk * nf_k;
    /*fprintf(stderr, "nf_jpkxnf_j=%f, nf_jpkxnf_k=%f\n", nf_jpkxnf_j, nf_jpkxnf_k);*/
    if ((nf_jpkxnf_j < GMX_FLOAT_EPS) || (nf_jpkxnf_k < GMX_FLOAT_EPS))
      return;
    cos_a = iprod(f_jpk, f_mj) / nf_jpkxnf_j;
    sin_a = norm(f_jpk_c_f_j) / nf_jpkxnf_j;
    cos_b = iprod(f_jpk, f_mk) / nf_jpkxnf_k;
    sin_b = norm(f_jpk_c_f_k) / nf_jpkxnf_k;
    /* in a triangle, known: length of one side and 2 angles; unknown: lengths of the 2 other sides */
    sinacosbpsinbcosa = sin_a * cos_b + sin_b * cos_a;
    /*fprintf(stderr, "sin_a=%f, cos_a=%f, sin_b=%f, cos_b=%f, sinacosbpsinbcosa=%f\n", sin_a, cos_a, sin_b, cos_b, sinacosbpsinbcosa);*/
    if (sinacosbpsinbcosa < GMX_FLOAT_EPS)
      return;
    nf_jpk_i = norm(f_jpk_i);
    nf_jpk_l = norm(f_jpk_l);
    nf_j_i = nf_jpk_i * sin_b / sinacosbpsinbcosa;
    nf_k_i = nf_jpk_i * sin_a / sinacosbpsinbcosa;
    nf_j_l = nf_jpk_l * sin_b / sinacosbpsinbcosa;
    nf_k_l = nf_jpk_l * sin_a / sinacosbpsinbcosa;
    /*fprintf(stderr, "nf_j_i=%f, nf_j_l=%f, nf_k_i=%f, nf_k_l=%f\n", nf_j_i, nf_j_l, nf_k_i, nf_k_l);*/
    //fprintf(stderr, "nf_k=%8f, nf_k_i+nf_k_l=%8f\n", nf_k, nf_k_i+nf_k_l);
    //fprintf(stderr, "nf_j=%8f, nf_j_i+nf_j_l=%8f\n", nf_j, nf_j_i+nf_j_l);
    /* make vectors from lengths */
    /* f_j_i and f_j_l are in the direction of f_j, f_k_i and f_k_l are in the direction of f_k */
    svmul(nf_j_i, uf_j, f_j_i);
    svmul(nf_j_l, uf_j, f_j_l);
    svmul(nf_k_i, uf_k, f_k_i);
    svmul(nf_k_l, uf_k, f_k_l);
    /* get f_j_k from difference */
    rvec_add(f_j_i, f_j_l, f_j_ipl);
    rvec_sub(f_mj, f_j_ipl, f_j_k);
    /* do the same for f_k_j, just to check by comparing with f_j_k */
    /*rvec_add(f_k_i, f_k_l, f_k_ipl);
    rvec_sub(f_mk, f_k_ipl, f_k_j);*/
    /* f_i is minus (f_j_i + f_k_i + f_l_i) because these are forces from i on these atoms, in the opposite direction from f_i */
    rvec_add(f_i, f_jpk_i, f_l_i);
    rvec_opp(f_l_i);
    /* do the same for f_il, just to check by comparing with f_il */
    /*rvec_add(f_l, f_jpk_l, f_i_l);
    rvec_opp(f_i_l);*/
    //fprintf(stderr, "f_i_l=%8f %8f %8f, f_l_i=%8f %8f %8f\n", f_i_l[0], f_i_l[1], f_i_l[2], f_l_i[0], f_l_i[1], f_l_i[2]);
    //fprintf(stderr, "nfi=%8f, nfij=%8f, nfik=%8f, nfil=%8f\n", norm(f_i), norm(f_ji), norm(f_ki), norm(f_il));
    /*forcemat_add_bonded(forcemat, j, k, f_jk, x, iBond);*/
    /*fprintf(debug, "f_ij=%8f, f_ik=%8f, f_il=%8f, f_kl=%8f, f_jl=%8f, f_jk=%8f\n", f_ij, f_ik, f_il, f_kl, f_jl, f_jk);*/
    /*rvec_sub(f_i, f_l, f_w);*/
    /*print_vec("f_w", f_w);*/
    /* !!! watch out, if uncommenting the below print_vec section, also uncomment the 2 sections above calculating f_k_j and f_i_l !!! */
    /*print_vec("f_i_l", f_i_l);
    print_vec("f_l_i", f_l_i);
    print_vec("f_j_i", f_j_i);
    print_vec("f_j_l", f_j_l);
    print_vec("f_k_i", f_k_i);
    print_vec("f_k_l", f_k_l);
    print_vec("f_jpk_i", f_jpk_i);
    print_vec("f_jpk_l", f_jpk_l);
    print_vec("f_j_k", f_j_k);
    print_vec("f_k_j", f_k_j);*/
    add_bonded(j, i, PF_INTER_DIHEDRAL, f_j_i);
    add_bonded(k, i, PF_INTER_DIHEDRAL, f_k_i);
    add_bonded(l, i, PF_INTER_DIHEDRAL, f_l_i);
    add_bonded(j, k, PF_INTER_DIHEDRAL, f_j_k);
    add_bonded(j, l, PF_INTER_DIHEDRAL, f_j_l);
    add_bonded(k, l, PF_INTER_DIHEDRAL, f_k_l);
}

void FDA::add_nonbonded(int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz)
{
  real pf_lj_residue, pf_coul_residue, pf_lj_coul;
  rvec pf_lj_atom_v, pf_lj_residue_v, pf_coul_atom_v, pf_coul_residue_v;

  /* first check that the interaction is interesting before doing expensive calculations and atom lookup*/
  if (!(fda_settings.type & PF_INTER_COULOMB))
    if (!(fda_settings.type & PF_INTER_LJ))
      return;
    else {
      add_nonbonded_single(i, j, PF_INTER_LJ, pf_lj, dx, dy, dz);
      return;
    }
  else
    if (!(fda_settings.type & PF_INTER_LJ)) {
      add_nonbonded_single(i, j, PF_INTER_COULOMB, pf_coul, dx, dy, dz);
      return;
    }

  if (!fda_settings.atoms_in_groups(i, j)) return;

  /* checking is symmetrical for atoms i and j; one of them has to be from g1, the other one from g2
   * however, if only residue_based_result_type is non-zero, atoms won't be initialized... so the conversion to residue numebers needs to be done here already;
   * the check below makes the atoms equivalent, make them always have the same order (i,j) and not (j,i) where i < j;
   * force is the force atom j exerts on atom i; if i and j are switched, the force goes in opposite direction
   * it's possible that i > j, but ri < rj, so the force has to be handled separately for each of them
   */
  if (residue_based.PF_or_PS_mode()) {
    /* the calling functions will not have i == j, but there is not such guarantee for ri and rj;
     * and it makes no sense to look at the interaction of a residue to itself
     */
    int ri = fda_settings.get_atom2residue(i);
    int rj = fda_settings.get_atom2residue(j);
    if (ri != rj) {
      if (ri > rj) {
        int tmp = rj; rj = ri; ri = tmp; // swap
        pf_lj_residue = -pf_lj;
        pf_coul_residue = -pf_coul;
      } else {
        pf_lj_residue = pf_lj;
        pf_coul_residue = pf_coul;
      }
      /* for detailed interactions, it's necessary to calculate and add separately, but for summed this can be simplified */
      switch(fda_settings.one_pair) {
        case OnePair::DETAILED:
          pf_coul_residue_v[0] = pf_coul_residue * dx;
          pf_coul_residue_v[1] = pf_coul_residue * dy;
          pf_coul_residue_v[2] = pf_coul_residue * dz;
          residue_based.distributed_forces.detailed[ri][rj][PF_INTER_COULOMB] += pf_coul_residue_v;
          pf_lj_residue_v[0] = pf_lj_residue * dx;
          pf_lj_residue_v[1] = pf_lj_residue * dy;
          pf_lj_residue_v[2] = pf_lj_residue * dz;
          residue_based.distributed_forces.detailed[ri][rj][PF_INTER_LJ] += pf_lj_residue_v;
          break;
        case OnePair::SUMMED:
          pf_lj_coul = pf_lj_residue + pf_coul_residue;
          pf_coul_residue_v[0] = pf_lj_coul * dx;
          pf_coul_residue_v[1] = pf_lj_coul * dy;
          pf_coul_residue_v[2] = pf_lj_coul * dz;
          residue_based.distributed_forces.summed[ri][rj] += pf_coul_residue_v;
          break;
        case OnePair::INVALID:
          break;
      }
    }
  }

  /* i & j as well as pf_lj & pf_coul are not used after this point, so it's safe to operate on their values directly */
  if (atom_based.PF_or_PS_mode()) {
    if (i > j) {
      int tmp = j; j = i; i = tmp; // swap
      pf_lj = -pf_lj;
      pf_coul = -pf_coul;
    }
    /* for detailed interactions, it's necessary to calculate and add separately, but for summed this can be simplified */
    switch(fda_settings.one_pair) {
      case OnePair::DETAILED:
        pf_coul_atom_v[0] = pf_coul * dx;
        pf_coul_atom_v[1] = pf_coul * dy;
        pf_coul_atom_v[2] = pf_coul * dz;
        atom_based.distributed_forces.detailed[i][j][PF_INTER_COULOMB] += pf_coul_atom_v;
        pf_lj_atom_v[0] = pf_lj * dx;
        pf_lj_atom_v[1] = pf_lj * dy;
        pf_lj_atom_v[2] = pf_lj * dz;
        atom_based.distributed_forces.detailed[i][j][PF_INTER_LJ] += pf_lj_atom_v;
        break;
      case OnePair::SUMMED:
        pf_lj_coul = pf_lj + pf_coul;
        pf_coul_atom_v[0] = pf_lj_coul * dx;
        pf_coul_atom_v[1] = pf_lj_coul * dy;
        pf_coul_atom_v[2] = pf_lj_coul * dz;
        atom_based.distributed_forces.summed[i][j] += pf_coul_atom_v;
        break;
      case OnePair::INVALID:
        break;
    }
  }
}

void FDA::add_virial(int ai, tensor v, real s)
{
	// Only symmetric tensor is used, therefore full multiplication is not as efficient
	//atom_vir[ai] += s * v;

    atom_vir[ai](XX, XX) += s * v[XX][XX];
    atom_vir[ai](YY, YY) += s * v[YY][YY];
    atom_vir[ai](ZZ, ZZ) += s * v[ZZ][ZZ];
    atom_vir[ai](XX, YY) += s * v[XX][YY];
    atom_vir[ai](XX, ZZ) += s * v[XX][ZZ];
    atom_vir[ai](YY, ZZ) += s * v[YY][ZZ];
}

void FDA::add_virial_bond(int ai, int aj, real f, real dx, real dy, real dz)
{
	tensor v;
	v[XX][XX] = dx * dx * f;
	v[YY][YY] = dy * dy * f;
	v[ZZ][ZZ] = dz * dz * f;
	v[XX][YY] = dx * dy * f;
	v[XX][ZZ] = dx * dz * f;
	v[YY][ZZ] = dy * dz * f;
	add_virial(ai, v, HALF);
	add_virial(aj, v, HALF);
}

void FDA::add_virial_angle(int ai, int aj, int ak,
  rvec r_ij, rvec r_kj, rvec f_i, rvec f_k)
{
  tensor v;
  v[XX][XX] = r_ij[XX] * f_i[XX] + r_kj[XX] * f_k[XX];
  v[YY][YY] = r_ij[YY] * f_i[YY] + r_kj[YY] * f_k[YY];
  v[ZZ][ZZ] = r_ij[ZZ] * f_i[ZZ] + r_kj[ZZ] * f_k[ZZ];
  v[XX][YY] = r_ij[XX] * f_i[YY] + r_kj[XX] * f_k[YY];
  v[XX][ZZ] = r_ij[XX] * f_i[ZZ] + r_kj[XX] * f_k[ZZ];
  v[YY][ZZ] = r_ij[YY] * f_i[ZZ] + r_kj[YY] * f_k[ZZ];
  add_virial(ai, v, THIRD);
  add_virial(aj, v, THIRD);
  add_virial(ak, v, THIRD);
}

void FDA::add_virial_dihedral(int i, int j, int k, int l,
  rvec f_i, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl)
{
  rvec r_lj;
  tensor v;
  rvec_sub(r_kj, r_kl, r_lj);
  v[XX][XX] = r_ij[XX] * f_i[XX] + r_kj[XX] * f_k[XX] + r_lj[XX] * f_l[XX];
  v[YY][YY] = r_ij[YY] * f_i[YY] + r_kj[YY] * f_k[YY] + r_lj[YY] * f_l[YY];
  v[ZZ][ZZ] = r_ij[ZZ] * f_i[ZZ] + r_kj[ZZ] * f_k[ZZ] + r_lj[ZZ] * f_l[ZZ];
  v[XX][YY] = r_ij[XX] * f_i[YY] + r_kj[XX] * f_k[YY] + r_lj[XX] * f_l[YY];
  v[XX][ZZ] = r_ij[XX] * f_i[ZZ] + r_kj[XX] * f_k[ZZ] + r_lj[XX] * f_l[ZZ];
  v[YY][ZZ] = r_ij[YY] * f_i[ZZ] + r_kj[YY] * f_k[ZZ] + r_lj[YY] * f_l[ZZ];
  add_virial(i, v, QUARTER);
  add_virial(j, v, QUARTER);
  add_virial(k, v, QUARTER);
  add_virial(l, v, QUARTER);
}

void FDA::write_frame(const rvec *x, gmx_mtop_t *mtop)
{
  rvec *com;
  atom_based.write_frame(x, nsteps);
  residue_based.write_frame(com, nsteps);
  ++nsteps;
}
