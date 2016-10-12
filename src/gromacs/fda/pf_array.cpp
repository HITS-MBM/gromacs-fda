/*
 * Handling of arrays of pairwise forces.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "pf_array.h"
#include "pf_array_detailed.h"
#include "pf_array_summed.h"
#include "pf_interactions.h"
#include "types/pf_array.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

static const real THIRD   =  1.0 / 3.0;
static const real QUARTER =  0.25;

/* check whether i and j are both in groups */
static inline gmx_bool pf_atoms_in_groups(int i, int j, t_pf_global *pf_global) {
  return ((pf_global->sys_in_g1[i] && pf_global->sys_in_g2[j]) || (pf_global->sys_in_g1[j] && pf_global->sys_in_g2[i]));
}

void pf_atom_add_bonded_nocheck(t_pf_global *pf_global, int i, int j, int type, rvec force) {
  t_pf_atoms *atoms;
  t_pf_atoms *residues;
  int ri = 0, rj = 0;       /* initialized to get rid of compiler warning, as they are only initialized later if ResidueBased is non-zero */
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

void pf_atom_add_bonded(t_pf_global *pf_global, int i, int j, int type, rvec force) {
  t_pf_atoms *atoms;
  t_pf_atoms *residues;
  int ai = 0, aj = 0;       /* initialized to get rid of compiler warning, as they are only initialized later if AtomBased is non-zero */
  int ri = 0, rj = 0;       /* initialized to get rid of compiler warning, as they are only initialized later if ResidueBased is non-zero */
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

void pf_atom_virial_angle(t_pf_global *pf_global, int ai, int aj, int ak,
	rvec r_ij, rvec r_kj, rvec f_i, rvec f_k)
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

void pf_atom_virial_dihedral(t_pf_global *pf_global, int i, int j, int k, int l,
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
	pf_atom_virial_add(pf_global, i, v, QUARTER);
	pf_atom_virial_add(pf_global, j, v, QUARTER);
	pf_atom_virial_add(pf_global, k, v, QUARTER);
	pf_atom_virial_add(pf_global, l, v, QUARTER);
}
