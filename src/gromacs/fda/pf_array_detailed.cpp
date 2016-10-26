/*
 * Functions related to detailed interactions between atoms.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#include <stddef.h>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "pf_array_detailed.h"
#include "pf_interactions.h"
#include "types/pf_array_detailed.h"
#include "fda.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

/* test for allocated to avoid creating a special function to just free arrays
 *
 * 2014-10-06 B.Doser @@ The array must not be allocated from scratch every md frame.
 * It's enough to zero a->len.
 */
static inline void pf_interaction_array_detailed_init(t_pf_interaction_array_detailed *a) {
  //if (a->allocated) sfree(a->array);
  //a->array = NULL;
  //a->allocated = 0;
  a->len = 0;
}

void pf_interactions_detailed_init(t_pf_interactions_detailed *interactions) {
  pf_interaction_array_detailed_init(&interactions->coulomb);
  pf_interaction_array_detailed_init(&interactions->lj);
  pf_interaction_array_detailed_init(&interactions->nb14);
  pf_interaction_array_detailed_init(&interactions->bonds);
  pf_interaction_array_detailed_init(&interactions->angles);
  pf_interaction_array_detailed_init(&interactions->dihedrals);
  pf_interaction_array_detailed_init(&interactions->polar);
}

void pf_atom_detailed_init(t_pf_atom_detailed *atom) {
  pf_interactions_detailed_init(&atom->interactions);
}

/* returns the t_pf_interaction_array_detailed corresponding to a certain interaction type */ 
t_pf_interaction_array_detailed *pf_interaction_array_by_type(t_pf_interactions_detailed *interactions, int type) {
  switch(type) {
    case PF_INTER_BOND:
      return &interactions->bonds;
    case PF_INTER_ANGLE:
      return &interactions->angles;
    case PF_INTER_DIHEDRAL:
      return &interactions->dihedrals;
    case PF_INTER_POLAR:
      return &interactions->polar;
    case PF_INTER_COULOMB:
      return &interactions->coulomb;
    case PF_INTER_LJ:
      return &interactions->lj;
    case PF_INTER_NB14:
      return &interactions->nb14;
    default:
      gmx_fatal(FARGS, "Unknown interaction type for pf: 0x%x.\n", type);
      /* not reached, but to avoid a compiler warning */
      return NULL;
  }
}

/* looks up an interaction with atom/residue jjnr and returns a pointer if found, otherwise NULL */
t_pf_interaction_detailed *pf_lookup_interaction_detailed(t_pf_interaction_array_detailed *ia, int jjnr) {
  int i;

  /* there is no check whether the array is allocated; if len is positive,
   * the array should have been allocated already
   */ 
  for (i = 0; (i < ia->len); i++)
    if (ia->array[i].jjnr == jjnr)
      return &ia->array[i];
  return NULL;
}

/* appends an interaction to the end of an interaction array */
void pf_interaction_array_detailed_append(t_pf_interaction_array_detailed *ia, int jjnr, rvec force) {
  /*fprintf(stderr, "pf_interaction_array_detailed_append ia=%p array=%p len=%d allocated=%d jjnr=%d force=%f,%f,%f\n", ia, ia->array, ia->len, ia->allocated, jjnr, force[0], force[1], force[2]);*/
  if (ia->len >= ia->allocated) {
    if (ia->array == NULL) {
      ia->allocated = PF_ARRAY_INITSIZE;
      snew(ia->array, ia->allocated);
    } else {
      ia->allocated = ia->allocated * OVER_ALLOC_FAC;
      srenew(ia->array, ia->allocated);
    }
  }
  ia->array[ia->len].jjnr = jjnr;
  copy_rvec(force, ia->array[ia->len].force);
  ia->len++;
}

/* looks for the interaction in the array corresponding to type; if found,
 * the force is added to the existing value; if not, a new interaction is
 * appended to the array
 */
void pf_atom_detailed_add(t_pf_atom_detailed *atom, int jjnr, int type, rvec force) {
  t_pf_interaction_array_detailed *ia;
  t_pf_interaction_detailed *p;

  /*fprintf(stderr, "pf_atom_add: jjnr=%d, type=0x%x, force=%f,%f,%f\n", jjnr, type, force[0], force[1], force[2]);*/
  ia = pf_interaction_array_by_type(&atom->interactions, type);
  p = pf_lookup_interaction_detailed(ia, jjnr);
  if (p == NULL) {
    pf_interaction_array_detailed_append(ia, jjnr, force);
  } else {
    rvec_inc(p->force, force);
  }
}

/* merge interactions of src into dst */
void pf_atom_detailed_merge(t_pf_atom_detailed *dst, t_pf_atom_detailed *src) {
  t_pf_interaction_array_detailed *ia_src, *ia_dst;
  t_pf_interaction_detailed *i_src, *i_dst;
  int type;
  int j;

  if (src->nr != dst->nr)
    gmx_fatal(FARGS, "Mismatch of PF atom detailed number: %d vs %d\n", src->nr, dst->nr);
  /* won't use pf_atom_detailed_add because the copying is done by interaction type, so the lookup can be done only once */
  for (type = 1; type < PF_INTER_ALL; type <<= 1) {
    ia_src = pf_interaction_array_by_type(&(src->interactions), type);
    ia_dst = pf_interaction_array_by_type(&(dst->interactions), type);
    for (j = 0; j < ia_src->len; j++){
      i_src = &(ia_src->array[j]);
      i_dst = pf_lookup_interaction_detailed(ia_dst, i_src->jjnr);
      if (i_dst == NULL)
        pf_interaction_array_detailed_append(ia_dst, i_src->jjnr, i_src->force);
      else
        rvec_inc(i_dst->force, i_src->force);
    }
  }
}

/* divides all forces by the given number; useful for calculating averages together with the pf_atom_*_merge() functions */
void pf_atom_detailed_real_divide(t_pf_atom_detailed *atom, real divisor) {
  t_pf_interaction_array_detailed *ia;
  t_pf_interaction_detailed *i;
  int type;
  int j;

  for (type = 1; type < PF_INTER_ALL; type <<= 1) {
    ia = pf_interaction_array_by_type(&(atom->interactions), type);
    for (j = 0; j < ia->len; j++)
      svdiv(divisor, ia->array[j].force);
  }
}
