/*
 * Handling of arrays of pairwise forces.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#include "fda.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "pf_array_detailed.h"
#include "pf_array_summed.h"
#include "pf_interactions.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

static const real HALF = 0.5;

/* check whether i and j are both in groups */
static inline gmx_bool pf_atoms_in_groups(int i, int j, struct t_pf_global *pf_global) {
  return ((pf_global->sys_in_g1[i] && pf_global->sys_in_g2[j]) || (pf_global->sys_in_g1[j] && pf_global->sys_in_g2[i]));
}

void pf_atom_add_nonbonded_single(struct t_pf_global *pf_global, int i, int j, int type, real force, real dx, real dy, real dz) {
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

void pf_atom_add_nonbonded(struct t_pf_global *pf_global, int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz) {
  t_pf_atoms *atoms;
  t_pf_atoms *residues;
  int ri = 0, rj = 0;
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
        int tmp = rj; rj = ri; ri = tmp; // swap
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
      int tmp = j; j = i; i = tmp; // swap
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

void pf_atom_virial_add(struct t_pf_global *pf_global, int ai, tensor v, real s)
{
    tensor *atom_vir = pf_global->atom_vir;
    atom_vir[ai][XX][XX] += s * v[XX][XX];
    atom_vir[ai][YY][YY] += s * v[YY][YY];
    atom_vir[ai][ZZ][ZZ] += s * v[ZZ][ZZ];
    atom_vir[ai][XX][YY] += s * v[XX][YY];
    atom_vir[ai][XX][ZZ] += s * v[XX][ZZ];
    atom_vir[ai][YY][ZZ] += s * v[YY][ZZ];
}

void pf_atom_virial_bond(struct t_pf_global *pf_global, int ai, int aj, real f, real dx, real dy, real dz)
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
