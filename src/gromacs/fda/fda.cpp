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
#include "pf_array_detailed.h"
#include "pf_array_summed.h"
#include "pf_interactions.h"
#include "pf_utils.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

static const real HALF    = 1.0 / 2.0;
static const real THIRD   = 1.0 / 3.0;
static const real QUARTER = 0.25;

FDA::FDA(FDASettings const& fda_settings)
 : fda_settings(fda_settings),
   atom_based_forces(ForceType::ATOMS, fda_settings),
   residue_based_forces(ForceType::RESIDUES, fda_settings),
   PFPS(0),
   VS(0),
   syslen_residues(0),
   sys_in_g1(nullptr),
   sys_in_g2(nullptr),
   atom2residue(nullptr),
   type(0),
   ofn_atoms(nullptr),
   of_atoms(nullptr),
   ofn_residues(nullptr),
   of_residues(nullptr),
   groupname(nullptr),
   nsteps_atoms(0),
   nsteps_residues(0),
   time_averages(nullptr),
   per_atom_real(nullptr),
   per_atom_real_int(nullptr),
   per_residue_real(nullptr),
   per_residue_real_int(nullptr),
   atom_vir(nullptr)
{
  /* check that there is an index file */
  if (!opt2bSet("-pfn", nfile, fnm))
	gmx_fatal(FARGS, "No index file (-pfn) for pairwise forces.\n");
  /* and finally read the index file and initialize the atoms lists */
  fprintf(stderr, "Pairwise forces for groups: %s and %s\n", g1name, g2name);

  PFPS = PF_or_PS_mode(fda_settings.atom_based_result_type) || PF_or_PS_mode(fda_settings.residue_based_result_type) ? 1 : 0;
  VS = VS_mode(fda_settings.atom_based_result_type) || VS_mode(fda_settings.residue_based_result_type) ? 1 : 0;

  /* if using compatibility mode, there should be only one group */
  if (compatibility_mode(fda_settings.atom_based_result_type) or compatibility_mode(fda_settings.residue_based_result_type)) {
	if (gmx_strcasecmp(g1name, g2name))
	  gmx_fatal(FARGS, "When using compat mode, the two group names should the identical.\n");
	else
	  groupname = gmx_strdup(g1name);
  }

  if ((stress_mode(fda_settings.atom_based_result_type) || stress_mode(fda_settings.residue_based_result_type)) && (OnePair != PF_ONEPAIR_SUMMED))
	gmx_fatal(FARGS, "Per atom data can only be computed from summed interactions.\n");

  pf_fill_atom2residue(this, top_global);	/* also fills syslen_residues */

  // allocates and fills with -1 the indexing table real atom nr. to pf number
  if (fda_settings.atom_based_result_type)
	atom_based_forces.forces.resize(fda_settings.syslen_atoms);
  if (fda_settings.residue_based_result_type)
	residue_based_forces.forces.resize(syslen_residues);

  if (PF_or_PS_mode(fda_settings.atom_based_result_type)) {
	pf_check_sys_in_g(this);
	pf_atoms_alloc(OnePair, atom_based_forces, fda_settings.syslen_atoms, "atoms");
	if (fda_settings.atom_based_result_type == ResultType::PUNCTUAL_STRESS) {
		pf_per_atom_real_init(&(per_atom_real), fda_settings.syslen_atoms, 0.0);
		ofn_atoms = gmx_strdup(opt2fn("-psa", nfile, fnm));
	} else {
		ofn_atoms = gmx_strdup(opt2fn("-pfa", nfile, fnm));
	}
	if (!ofn_atoms)
	  gmx_fatal(FARGS, "No file for writing out atom-based values.\n");
  }

  if (PF_or_PS_mode(fda_settings.residue_based_result_type)) {
	pf_check_sys_in_g(this);
	pf_atoms_alloc(OnePair, residue_based_forces, syslen_residues, "residues");
	if (fda_settings.residue_based_result_type == ResultType::PUNCTUAL_STRESS) {
		pf_per_atom_real_init(&(per_residue_real), syslen_residues, 0.0);
		ofn_residues = gmx_strdup(opt2fn("-psr", nfile, fnm));
	} else {
		ofn_residues = gmx_strdup(opt2fn("-pfr", nfile, fnm));
	}
	if (!ofn_residues)
	  gmx_fatal(FARGS, "No file for writing out residue-based values.\n");
  }

  if (VS_mode(fda_settings.atom_based_result_type)) {
	snew(atom_vir, fda_settings.syslen_atoms);
	if (fda_settings.atom_based_result_type == ResultType::VIRIAL_STRESS)
		ofn_atoms = gmx_strdup(opt2fn("-vsa", nfile, fnm));
	if (fda_settings.atom_based_result_type == ResultType::VIRIAL_STRESS_VON_MISES)
		ofn_atoms = gmx_strdup(opt2fn("-vma", nfile, fnm));
	if (!ofn_atoms) gmx_fatal(FARGS, "No file for writing out virial stress.\n");
  }

  /* initialization of time averages */
  snew(time_averages, 1);
  ITYPE("time_averages_period", time_averages->period, 1);
  if (time_averages->period < 0)
	gmx_fatal(FARGS, "Invalid value for time_averages_period: %d\n", time_averages->period);
  time_averages->steps = 0;
  if (time_averages->period != 1) {
	if (OnePair != PF_ONEPAIR_SUMMED)
	  gmx_fatal(FARGS, "Can only save scalar time averages from summed interactions.\n");
	/* for atoms/residues, pf_atoms_init is called for each step from md.c, but for time averages the initialization needs to be done only at the begining and after every data writing */
	if (PF_or_PS_mode(fda_settings.atom_based_result_type)) {
	  if (!((compatibility_mode(fda_settings.atom_based_result_type)) || (fda_settings.atom_based_result_type == ResultType::PAIRWISE_FORCES_SCALAR)))
		gmx_fatal(FARGS, "Can only use time averages with scalar or compatibility output.\n");
	  pf_atoms_scalar_alloc(atom_based_forces, fda_settings.syslen_atoms, "atoms - time averages");
	  pf_atoms_scalar_init(atom_based_forces);
	}
	if (PF_or_PS_mode(fda_settings.residue_based_result_type)) {
	  if (!((compatibility_mode(fda_settings.residue_based_result_type)) || (fda_settings.residue_based_result_type == ResultType::PAIRWISE_FORCES_SCALAR)))
		gmx_fatal(FARGS, "Can only use time averages with scalar or compatibility output.\n");
	  pf_atoms_scalar_alloc(residue_based_forces, syslen_residues, "residues - time averages");
	  pf_atoms_scalar_init(residue_based_forces);
	  snew(time_averages->com, syslen_residues);
	  clear_rvecs(syslen_residues, time_averages->com);
	}
  }
}

void FDA::add_bonded_nocheck(int i, int j, int type, rvec force)
{
  /* checking is symmetrical for atoms i and j; one of them has to be from g1, the other one from g2;
   * the check below makes the atoms equivalent, make them always have the same order (i,j) and not (j,i) where i < j;
   * force is the force atom j exerts on atom i; if i and j are switched, the force goes in opposite direction
   * it's possible that i > j, but ri < rj, so the force has to be handled separately for each of them
   * 
   * the logic is a bit complicated by the fact that atom_based_result_type and residue_based_result_type are independent;
   * if residue_based_result_type part is done first, the atom_based_result_type part can use i/j/force directly, without saving them
   * first in intermediate variables, as the initial values of i/j/force are no longer needed; if atom_based_result_type
   * is done first (the original code), i/j/force are needed for the later atom->residue mapping
   * and saving in intermediate variables is needed
   */
  if (residue_based_forces.PF_or_PS_mode()) {
    /* the calling functions will not have i == j, but there is not such guarantee for ri and rj;
     * and it makes no sense to look at the interaction of a residue to itself
     */
    int ri = atom2residue[i];
    int rj = atom2residue[j];
    //fprintf(stderr, "pf_atom_add_bonded_nocheck: i=%d, j=%d, ri=%d, rj=%d, type=%d\n", i, j, ri, rj, type);
    rvec force_residue;
    if (ri != rj) {
      switch(fda_settings.one_pair) {
        case OnePair::DETAILED:
          if (ri > rj) {
            int_swap(&ri, &rj);
            clear_rvec(force_residue);
            rvec_dec(force_residue, force);
            pf_atom_detailed_add(&residue_based_forces->detailed[residue_based_forces->sys2pf[ri]], rj, type, force_residue);
          } else {
            pf_atom_detailed_add(&residue_based_forces->detailed[residue_based_forces->sys2pf[ri]], rj, type, force);
          }
          break;
        case OnePair::SUMMED:
          if (ri > rj) {
            int_swap(&ri, &rj);
            clear_rvec(force_residue);
            rvec_dec(force_residue, force);
            pf_atom_summed_add(&residue_based_forces->summed[residue_based_forces->sys2pf[ri]], rj, type, force_residue);
          } else {
            pf_atom_summed_add(&residue_based_forces->summed[residue_based_forces->sys2pf[ri]], rj, type, force);
          }
          break;
        default:
          break;
      }
    }
  }

  if (atom_based_forces.PF_or_PS_mode()) {
    //fprintf(stderr, "pf_atom_add_bonded_nocheck: i=%d, j=%d, type=%d\n", i, j, type);
    if (i > j) {
      int_swap(&i, &j);
      rvec_opp(force);
    }
    switch(fda_settings.one_pair) {
      case OnePair::DETAILED:
        pf_atom_detailed_add(&atom_based_forces->detailed[atom_based_forces->sys2pf[i]], j, type, force);
        break;
      case OnePair::SUMMED:
        pf_atom_summed_add(&atom_based_forces->summed[atom_based_forces->sys2pf[i]], j, type, force);
        break;
      default:
        break;
    }
  }
}

void FDA::add_bonded(int i, int j, int type, rvec force)
{
  // leave early if the interaction is not interesting
  if (!(this->type & type)) return;
  if (!this->atoms_in_groups(i, j)) {
	//fprintf(stderr, "Warning: Unneeded pair in pf_atom_add_bonded i=%i, j=%i\n", i, j); fflush(stderr);
    return;
  }
  add_bonded_nocheck(i, j, type, force);
}

void FDA::add_nonbonded_single(int i, int j, int type, real force, real dx, real dy, real dz)
{
  rvec force_v;			/* vector force for interaction */

  /* first check that the interaction is interesting before doing computation and lookups */
  if (!(type & type))
    return;
  if (!atoms_in_groups(i, j)) {
	//fprintf(stderr, "Warning: Unneeded pair in pf_atom_add_nonbonded_single i=%i, j=%i\n", i, j); fflush(stderr);
    return;
  }

  force_v[0] = force * dx;
  force_v[1] = force * dy;
  force_v[2] = force * dz;
  add_bonded_nocheck(i, j, type, force_v);
}

void FDA::add_nonbonded(int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz)
{
  int ri = 0, rj = 0;
  real pf_lj_residue, pf_coul_residue, pf_lj_coul;
  rvec pf_lj_atom_v, pf_lj_residue_v, pf_coul_atom_v, pf_coul_residue_v;

  /* first check that the interaction is interesting before doing expensive calculations and atom lookup*/
  if (!(type & PF_INTER_COULOMB))
    if (!(type & PF_INTER_LJ))
      return;
    else {
      add_nonbonded_single(i, j, PF_INTER_LJ, pf_lj, dx, dy, dz);
      return;
    }
  else
    if (!(type & PF_INTER_LJ)) {
      add_nonbonded_single(i, j, PF_INTER_COULOMB, pf_coul, dx, dy, dz);
      return;
    }
  if (!atoms_in_groups(i, j)) {
	//fprintf(stderr, "Warning: Unneeded pair in pf_atom_add_nonbonded i=%i, j=%i\n", i, j); fflush(stderr);
    return;
  }

  /*fprintf(stderr, "Nonbonded interaction: ii=%d, jnr=%d\n", ii, jnr);*/

  /* checking is symmetrical for atoms i and j; one of them has to be from g1, the other one from g2
   * however, if only residue_based_result_type is non-zero, atoms won't be initialized... so the conversion to residue numebers needs to be done here already;
   * the check below makes the atoms equivalent, make them always have the same order (i,j) and not (j,i) where i < j;
   * force is the force atom j exerts on atom i; if i and j are switched, the force goes in opposite direction
   * it's possible that i > j, but ri < rj, so the force has to be handled separately for each of them
   */
  if (pf_file_out_PF_or_PS(fda_settings.residue_based_result_type)) {
    /* the calling functions will not have i == j, but there is not such guarantee for ri and rj;
     * and it makes no sense to look at the interaction of a residue to itself
     */
    ri = atom2residue[i];
    rj = atom2residue[j];
    //fprintf(stderr, "pf_atom_add_nonbonded: i=%d, j=%d, ri=%d, rj=%d\n", i, j, ri, rj);
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
      switch(OnePair) {
        case PF_ONEPAIR_DETAILED:
          pf_coul_residue_v[0] = pf_coul_residue * dx;
          pf_coul_residue_v[1] = pf_coul_residue * dy;
          pf_coul_residue_v[2] = pf_coul_residue * dz;
          pf_lj_residue_v[0] = pf_lj_residue * dx;
          pf_lj_residue_v[1] = pf_lj_residue * dy;
          pf_lj_residue_v[2] = pf_lj_residue * dz;
          pf_atom_detailed_add(&residue_based_forces->detailed[residue_based_forces->sys2pf[ri]], rj, PF_INTER_LJ, pf_lj_residue_v);
          pf_atom_detailed_add(&residue_based_forces->detailed[residue_based_forces->sys2pf[ri]], rj, PF_INTER_COULOMB, pf_coul_residue_v);
          break;
        case PF_ONEPAIR_SUMMED:
          pf_lj_coul = pf_lj_residue + pf_coul_residue;
          pf_coul_residue_v[0] = pf_lj_coul * dx;
          pf_coul_residue_v[1] = pf_lj_coul * dy;
          pf_coul_residue_v[2] = pf_lj_coul * dz;
          pf_atom_summed_add(&residue_based_forces->summed[residue_based_forces->sys2pf[ri]], rj, PF_INTER_LJ | PF_INTER_COULOMB, pf_coul_residue_v);
          break;
        default:
          break;
      }
    }
  }

  /* i & j as well as pf_lj & pf_coul are not used after this point, so it's safe to operate on their values directly */
  if (pf_file_out_PF_or_PS(settings.atom_based_result_type)) {
    //fprintf(stderr, "pf_atom_add_nonbonded: i=%d, j=%d\n", i, j);
    if (i > j) {
      int tmp = j; j = i; i = tmp; // swap
      pf_lj = -pf_lj;
      pf_coul = -pf_coul;
    }
    /* for detailed interactions, it's necessary to calculate and add separately, but for summed this can be simplified */
    switch(OnePair) {
      case PF_ONEPAIR_DETAILED:
        pf_coul_atom_v[0] = pf_coul * dx;
        pf_coul_atom_v[1] = pf_coul * dy;
        pf_coul_atom_v[2] = pf_coul * dz;
        pf_lj_atom_v[0] = pf_lj * dx;
        pf_lj_atom_v[1] = pf_lj * dy;
        pf_lj_atom_v[2] = pf_lj * dz;
        pf_atom_detailed_add(&atom_based_forces->detailed[atom_based_forces->sys2pf[i]], j, PF_INTER_LJ, pf_lj_atom_v);
        pf_atom_detailed_add(&atom_based_forces->detailed[atom_based_forces->sys2pf[i]], j, PF_INTER_COULOMB, pf_coul_atom_v);
        break;
      case PF_ONEPAIR_SUMMED:
        pf_lj_coul = pf_lj + pf_coul;
        pf_coul_atom_v[0] = pf_lj_coul * dx;
        pf_coul_atom_v[1] = pf_lj_coul * dy;
        pf_coul_atom_v[2] = pf_lj_coul * dz;
        pf_atom_summed_add(&atom_based_forces->summed[atom_based_forces->sys2pf[i]], j, PF_INTER_LJ | PF_INTER_COULOMB, pf_coul_atom_v);
        break;
      default:
        break;
    }
  }
}

void FDA::add_virial(int ai, tensor v, real s)
{
    atom_vir[ai][XX][XX] += s * v[XX][XX];
    atom_vir[ai][YY][YY] += s * v[YY][YY];
    atom_vir[ai][ZZ][ZZ] += s * v[ZZ][ZZ];
    atom_vir[ai][XX][YY] += s * v[XX][YY];
    atom_vir[ai][XX][ZZ] += s * v[XX][ZZ];
    atom_vir[ai][YY][ZZ] += s * v[YY][ZZ];
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
