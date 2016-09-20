/*
 * Handling of arrays of pairwise forces.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#include "gromacs/legacyheaders/pf_array_detailed.h"
#include "gromacs/legacyheaders/pf_array_summed.h"
#include "gromacs/legacyheaders/pf_interactions.h"
#include "gromacs/legacyheaders/types/pf_array.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

static const real HALF    =  0.5;
static const real THIRD   =  1.0 / 3.0;
static const real QUARTER =  0.25;

/* check whether i and j are both in groups */
static inline gmx_bool pf_atoms_in_groups(int i, int j, t_pf_global *pf_global) {
  return ((pf_global->sys_in_g1[i] && pf_global->sys_in_g2[j]) || (pf_global->sys_in_g1[j] && pf_global->sys_in_g2[i]));
}

static inline void pf_atom_add_bonded_nocheck(t_pf_global *pf_global, atom_id i, atom_id j, int type, rvec force) {
  t_pf_atoms *atoms;
  t_pf_atoms *residues;
  atom_id ri = 0, rj = 0;       /* initialized to get rid of compiler warning, as they are only initialized later if ResidueBased is non-zero */
  rvec force_residue;

  atoms = pf_global->atoms;
  residues = pf_global->residues;
  //fprintf(stderr, "pf_atom_add_bonded_nocheck: adding i=%d, j=%d, type=%d\n", i, j, type);

  /* checking is symmetrical for atoms i and j; one of them has to be from g1, the other one from g2;
   * the check below makes the atoms equivalent, make them always have the same order (i,j) and not (j,i) where i < j;
   * force is the force atom j exerts on atom i; if i and j are switched, the force goes in opposite direction
   * it's possible that i > j, but ri < rj, so the force has to be handled separately for each of them
   * 
   * the logic is a bit complicated by the fact that AtomBased and ResidueBased are independent;
   * if ResidueBased part is done first, the AtomBased part can use i/j/force directly, without saving them
   * first in intermediate variables, as the initial values of i/j/force are no longer needed; if AtomBased
   * is done first (the original code), i/j/force are needed for the later atom->residue mapping
   * and saving in intermediate variables is needed
   */
  if (pf_file_out_PF_or_PS(pf_global->ResidueBased)) {
    /* the calling functions will not have i == j, but there is not such guarantee for ri and rj;
     * and it makes no sense to look at the interaction of a residue to itself
     */
    ri = pf_global->atom2residue[i];
    rj = pf_global->atom2residue[j];
    //fprintf(stderr, "pf_atom_add_bonded_nocheck: i=%d, j=%d, ri=%d, rj=%d, type=%d\n", i, j, ri, rj, type);
    if (ri != rj) {
      switch(pf_global->OnePair) {
        case PF_ONEPAIR_DETAILED:
          if (ri > rj) {
            int_swap(&ri, &rj);
            clear_rvec(force_residue);
            rvec_dec(force_residue, force);
            pf_atom_detailed_add(&residues->detailed[residues->sys2pf[ri]], rj, type, force_residue);
          } else {
            pf_atom_detailed_add(&residues->detailed[residues->sys2pf[ri]], rj, type, force);
          }
          break;
        case PF_ONEPAIR_SUMMED:
          if (ri > rj) {
            int_swap(&ri, &rj);
            clear_rvec(force_residue);
            rvec_dec(force_residue, force);
            pf_atom_summed_add(&residues->summed[residues->sys2pf[ri]], rj, type, force_residue);
          } else {
            pf_atom_summed_add(&residues->summed[residues->sys2pf[ri]], rj, type, force);
          }
          break;
        default:
          break;
      }
    }
  }

  if (pf_file_out_PF_or_PS(pf_global->AtomBased)) {
    //fprintf(stderr, "pf_atom_add_bonded_nocheck: i=%d, j=%d, type=%d\n", i, j, type);
    if (i > j) {
      int_swap(&i, &j);
      rvec_opp(force);
    }
    switch(pf_global->OnePair) {
      case PF_ONEPAIR_DETAILED:
        pf_atom_detailed_add(&atoms->detailed[atoms->sys2pf[i]], j, type, force);
        break;
      case PF_ONEPAIR_SUMMED:
        pf_atom_summed_add(&atoms->summed[atoms->sys2pf[i]], j, type, force);
        break;
      default:
        break;
    }
  }
}

void pf_atom_add_bonded(t_pf_global *pf_global, atom_id i, atom_id j, int type, rvec force) {
  t_pf_atoms *atoms;
  t_pf_atoms *residues;
  atom_id ai = 0, aj = 0;       /* initialized to get rid of compiler warning, as they are only initialized later if AtomBased is non-zero */
  atom_id ri = 0, rj = 0;       /* initialized to get rid of compiler warning, as they are only initialized later if ResidueBased is non-zero */
  rvec force_atom, force_residue;
  gmx_bool atom_add = FALSE;
  gmx_bool residue_add = FALSE;

  /* leave early if the interaction is not interesting */
  if (!(pf_global->type & type))
    return;
  if (!pf_atoms_in_groups(i, j, pf_global)) {
	//fprintf(stderr, "Warning: Unneeded pair in pf_atom_add_bonded i=%i, j=%i\n", i, j); fflush(stderr);
    return;
  }

  pf_atom_add_bonded_nocheck(pf_global, i, j, type, force);
}

/* add a particular type of nonbonded interaction for the kernels where only one type of interaction is computed;
 * force is passed as scalar along with the distance vector (as dx, dy, dz) from which the vector force is
 * computed, the same way it's done in the nonbonded kernels
 */
void pf_atom_add_nonbonded_single(t_pf_global *pf_global, atom_id i, atom_id j, int type, real force, real dx, real dy, real dz) {
  rvec force_v;			/* vector force for interaction */

  /* first check that the interaction is interesting before doing computation and lookups */
  if (!(pf_global->type & type))
    return;
  if (!pf_atoms_in_groups(i, j, pf_global)) {
	//fprintf(stderr, "Warning: Unneeded pair in pf_atom_add_nonbonded_single i=%i, j=%i\n", i, j); fflush(stderr);
    return;
  }

  force_v[0] = force * dx;
  force_v[1] = force * dy;
  force_v[2] = force * dz;
  pf_atom_add_bonded_nocheck(pf_global, i, j, type, force_v);
}

/* add a nonbonded interaction for kernels where both Coulomb and LJ are computed;
 * this is more efficient than calling the previous one twice because some of the tests are made only once;
 * forces are passed as scalars along with the distance vector (as dx, dy, dz) from which the vector forces are
 * computed, the same way it's done in the nonbonded kernels
 */
void pf_atom_add_nonbonded(t_pf_global *pf_global, atom_id i, atom_id j, real pf_coul, real pf_lj, real dx, real dy, real dz) {
  t_pf_atoms *atoms;
  t_pf_atoms *residues;
  atom_id ri = 0, rj = 0;
  real pf_lj_residue, pf_coul_residue, pf_lj_coul;
  rvec pf_lj_atom_v, pf_lj_residue_v, pf_coul_atom_v, pf_coul_residue_v;

  /* first check that the interaction is interesting before doing expensive calculations and atom lookup*/
  if (!(pf_global->type & PF_INTER_COULOMB))
    if (!(pf_global->type & PF_INTER_LJ))
      return;
    else {
      pf_atom_add_nonbonded_single(pf_global, i, j, PF_INTER_LJ, pf_lj, dx, dy, dz);
      return;
    }
  else
    if (!(pf_global->type & PF_INTER_LJ)) {
      pf_atom_add_nonbonded_single(pf_global, i, j, PF_INTER_COULOMB, pf_coul, dx, dy, dz);
      return;
    }
  if (!pf_atoms_in_groups(i, j, pf_global)) {
	//fprintf(stderr, "Warning: Unneeded pair in pf_atom_add_nonbonded i=%i, j=%i\n", i, j); fflush(stderr);
    return;
  }

  atoms = pf_global->atoms;
  residues = pf_global->residues;
  /*fprintf(stderr, "Nonbonded interaction: ii=%d, jnr=%d\n", ii, jnr);*/

  /* checking is symmetrical for atoms i and j; one of them has to be from g1, the other one from g2
   * however, if only ResidueBased is non-zero, pf_global->atoms won't be initialized... so the conversion to residue numebers needs to be done here already;
   * the check below makes the atoms equivalent, make them always have the same order (i,j) and not (j,i) where i < j;
   * force is the force atom j exerts on atom i; if i and j are switched, the force goes in opposite direction
   * it's possible that i > j, but ri < rj, so the force has to be handled separately for each of them
   */
  if (pf_file_out_PF_or_PS(pf_global->ResidueBased)) {
    /* the calling functions will not have i == j, but there is not such guarantee for ri and rj;
     * and it makes no sense to look at the interaction of a residue to itself
     */
    ri = pf_global->atom2residue[i];
    rj = pf_global->atom2residue[j];
    //fprintf(stderr, "pf_atom_add_nonbonded: i=%d, j=%d, ri=%d, rj=%d\n", i, j, ri, rj);
    if (ri != rj) {
      if (ri > rj) {
        int_swap(&ri, &rj);
        pf_lj_residue = -pf_lj;
        pf_coul_residue = -pf_coul;
      } else {
        pf_lj_residue = pf_lj;
        pf_coul_residue = pf_coul;
      }
      /* for detailed interactions, it's necessary to calculate and add separately, but for summed this can be simplified */
      switch(pf_global->OnePair) {
        case PF_ONEPAIR_DETAILED:
          pf_coul_residue_v[0] = pf_coul_residue * dx;
          pf_coul_residue_v[1] = pf_coul_residue * dy;
          pf_coul_residue_v[2] = pf_coul_residue * dz;
          pf_lj_residue_v[0] = pf_lj_residue * dx;
          pf_lj_residue_v[1] = pf_lj_residue * dy;
          pf_lj_residue_v[2] = pf_lj_residue * dz;
          pf_atom_detailed_add(&residues->detailed[residues->sys2pf[ri]], rj, PF_INTER_LJ, pf_lj_residue_v);
          pf_atom_detailed_add(&residues->detailed[residues->sys2pf[ri]], rj, PF_INTER_COULOMB, pf_coul_residue_v);
          break;
        case PF_ONEPAIR_SUMMED:
          pf_lj_coul = pf_lj_residue + pf_coul_residue;
          pf_coul_residue_v[0] = pf_lj_coul * dx;
          pf_coul_residue_v[1] = pf_lj_coul * dy;
          pf_coul_residue_v[2] = pf_lj_coul * dz;
          pf_atom_summed_add(&residues->summed[residues->sys2pf[ri]], rj, PF_INTER_LJ | PF_INTER_COULOMB, pf_coul_residue_v);
          break;
        default:
          break;
      }
    }
  }

  /* i & j as well as pf_lj & pf_coul are not used after this point, so it's safe to operate on their values directly */
  if (pf_file_out_PF_or_PS(pf_global->AtomBased)) {
    //fprintf(stderr, "pf_atom_add_nonbonded: i=%d, j=%d\n", i, j);
    if (i > j) {
      int_swap(&i, &j);
      pf_lj = -pf_lj;
      pf_coul = -pf_coul;
    }
    /* for detailed interactions, it's necessary to calculate and add separately, but for summed this can be simplified */
    switch(pf_global->OnePair) {
      case PF_ONEPAIR_DETAILED:
        pf_coul_atom_v[0] = pf_coul * dx;
        pf_coul_atom_v[1] = pf_coul * dy;
        pf_coul_atom_v[2] = pf_coul * dz;
        pf_lj_atom_v[0] = pf_lj * dx;
        pf_lj_atom_v[1] = pf_lj * dy;
        pf_lj_atom_v[2] = pf_lj * dz;
        pf_atom_detailed_add(&atoms->detailed[atoms->sys2pf[i]], j, PF_INTER_LJ, pf_lj_atom_v);
        pf_atom_detailed_add(&atoms->detailed[atoms->sys2pf[i]], j, PF_INTER_COULOMB, pf_coul_atom_v);
        break;
      case PF_ONEPAIR_SUMMED:
        pf_lj_coul = pf_lj + pf_coul;
        pf_coul_atom_v[0] = pf_lj_coul * dx;
        pf_coul_atom_v[1] = pf_lj_coul * dy;
        pf_coul_atom_v[2] = pf_lj_coul * dz;
        pf_atom_summed_add(&atoms->summed[atoms->sys2pf[i]], j, PF_INTER_LJ | PF_INTER_COULOMB, pf_coul_atom_v);
        break;
      default:
        break;
    }
  }
}

/**
 * The atom virial can be expressed as a 6-real tensor, as it's symmetric.
 * To avoid defining a new tensor type, the 9-real tensor is used instead.
 */
void pf_atom_virial_add(t_pf_global *pf_global, int ai, tensor v, real s)
{
    tensor *atom_vir = pf_global->atom_vir;
    atom_vir[ai][XX][XX] += s * v[XX][XX];
    atom_vir[ai][YY][YY] += s * v[YY][YY];
    atom_vir[ai][ZZ][ZZ] += s * v[ZZ][ZZ];
    atom_vir[ai][XX][YY] += s * v[XX][YY];
    atom_vir[ai][XX][ZZ] += s * v[XX][ZZ];
    atom_vir[ai][YY][ZZ] += s * v[YY][ZZ];
}

/**
 * Origin on j, but for 2 atoms it doesn't matter.
 */
void pf_atom_virial_bond(t_pf_global *pf_global, int ai, int aj, real f, real dx, real dy, real dz)
{
	tensor v;
	v[XX][XX] = dx * dx * f;
	v[YY][YY] = dy * dy * f;
	v[ZZ][ZZ] = dz * dz * f;
	v[XX][YY] = dx * dy * f;
	v[XX][ZZ] = dx * dz * f;
	v[YY][ZZ] = dy * dz * f;
	pf_atom_virial_add(pf_global, ai, v, HALF);
	pf_atom_virial_add(pf_global, aj, v, HALF);
}

/**
 * Translate to origin on the middle (j) atom:
 * vir = ri*Fi + rj*Fj + rk*Fk
 *     = (ri-rj)*Fi + (rk-rj)*Fk
 *     = r_ij[dim1]*f_i[dim2] + r_kj[dim1]*f_k[dim2]
 */
void pf_atom_virial_angle(t_pf_global *pf_global, int ai, int aj, int ak, rvec r_ij, rvec r_kj, rvec f_i, rvec f_k)
{
	tensor v;
	v[XX][XX] = r_ij[XX] * f_i[XX] + r_kj[XX] * f_k[XX];
	v[YY][YY] = r_ij[YY] * f_i[YY] + r_kj[YY] * f_k[YY];
	v[ZZ][ZZ] = r_ij[ZZ] * f_i[ZZ] + r_kj[ZZ] * f_k[ZZ];
	v[XX][YY] = r_ij[XX] * f_i[YY] + r_kj[XX] * f_k[YY];
	v[XX][ZZ] = r_ij[XX] * f_i[ZZ] + r_kj[XX] * f_k[ZZ];
	v[YY][ZZ] = r_ij[YY] * f_i[ZZ] + r_kj[YY] * f_k[ZZ];
	pf_atom_virial_add(pf_global, ai, v, THIRD);
	pf_atom_virial_add(pf_global, aj, v, THIRD);
	pf_atom_virial_add(pf_global, ak, v, THIRD);
}

/**
 * Translate to origin on the second (j) atom:
 * vir = ri*Fi + rj*Fj + rk*Fk + rl*Fl
 *     = (ri-rj)*Fi + (rk-rj)*Fk + (rl-rj)*Fl
 *      = (ri-rj)*Fi + (rk-rj)*Fk + ((rl-rk) + (rk-rj))*Fl
 *      = r_ij[dim1]*f_i[dim2] + r_kj[dim1]*f_k[dim2] + (r_kj-r_kl)[dim1]*f_l[dim2]
 */
void pf_atom_virial_dihedral(t_pf_global *pf_global, int i, int j, int k, int l, rvec f_i, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl)
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
	pf_atom_virial_add(pf_global, i, v, QUARTER);
	pf_atom_virial_add(pf_global, j, v, QUARTER);
	pf_atom_virial_add(pf_global, k, v, QUARTER);
	pf_atom_virial_add(pf_global, l, v, QUARTER);
}
