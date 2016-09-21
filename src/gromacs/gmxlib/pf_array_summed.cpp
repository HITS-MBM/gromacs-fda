/*
 * Functions related to summed interactions between atoms.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#include <stddef.h>
#include "gromacs/legacyheaders/types/pf_array_summed.h"
#include "gromacs/legacyheaders/types/pf_array.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

/* test for allocated to avoid creating a special function to just free arrays
 *
 * 2014-10-06 B.Doser @@ The array must not be allocated from scratch every md frame.
 * It's enough to zero a->len.
 */
static inline void pf_interaction_array_summed_init(t_pf_interaction_array_summed *a) {
  //if (a->allocated) sfree(a->array);
  //a->array = NULL;
  //a->allocated = 0;
  a->len = 0;
}

void pf_atom_summed_init(t_pf_atom_summed *atom) {
  pf_interaction_array_summed_init(&atom->interactions);
}

t_pf_interaction_summed *pf_lookup_interaction_summed(t_pf_interaction_array_summed *ia, int jjnr) {
  atom_id i;

  /* there is no check whether the array is allocated; if len is positive,
   * the array should have been allocated already
   */ 
  for (i = 0; (i < ia->len); i++)
    if (ia->array[i].jjnr == jjnr)
      return &ia->array[i];
  return NULL;
}

/* appends an interaction to the end of an interaction array */
void pf_interaction_array_summed_append(t_pf_interaction_array_summed *ia, atom_id jjnr, int type, rvec force) {
  //fprintf(stderr, "pf_interaction_array_summed_append ia_summed=%p array=%p len=%d allocated=%d jjnr=%d force=%f,%f,%f\n", ia, ia->array, ia->len, ia->allocated, jjnr, force[0], force[1], force[2]);
  if (ia->len >= ia->allocated) {
    if (ia->array == NULL) {
  	  //fprintf(stderr, "pf_interaction_array_summed.array is allocated, old: %i, new: %i\n", ia->allocated, PF_ARRAY_INITSIZE);
      ia->allocated = PF_ARRAY_INITSIZE;
      snew(ia->array, ia->allocated);
    } else {
      //fprintf(stderr, "pf_interaction_array_summed.array is reallocated, old: %i, new: %i\n", ia->allocated, (int)(ia->allocated * OVER_ALLOC_FAC));
      ia->allocated = ia->allocated * OVER_ALLOC_FAC;
      srenew(ia->array, ia->allocated);
    }
  }
  ia->array[ia->len].jjnr = jjnr;
  ia->array[ia->len].type = type;
  copy_rvec(force, ia->array[ia->len].force);
  ia->len++;
}

/* looks for the interaction in the array corresponding to type; if found,
 * the force is added to the existing value; if not, a new interaction is
 * appended to the array
 */
void pf_atom_summed_add(t_pf_atom_summed *atom, atom_id jjnr, int type, rvec force) {
  t_pf_interaction_summed *p;

  //fprintf(stderr, "pf_atom_summed_add: jjnr=%d, type=%d, force=%e,%e,%e\n", jjnr, type, force[0], force[1], force[2]);
  p = pf_lookup_interaction_summed(&atom->interactions, jjnr);
  if (p == NULL) {
    pf_interaction_array_summed_append(&atom->interactions, jjnr, type, force);
  } else {
    rvec_inc(p->force, force);
    p->type |= type;
  }
}

/* merge interactions of src into dst */
void pf_atom_summed_merge(t_pf_atom_summed *dst, t_pf_atom_summed *src) {
  t_pf_interaction_array_summed *ia;
  t_pf_interaction_summed *i;
  atom_id j;

  if (src->nr != dst->nr)
    gmx_fatal(FARGS, "Mismatch of PF atom summed number: %d vs %d\n", src->nr, dst->nr);
  ia = &(src->interactions);
  for (j = 0; j < ia->len; j++) {
    i = &(ia->array[j]);
    pf_atom_summed_add(dst, i->jjnr, i->type, i->force);
  }
}

/* divides all forces by the given number; useful for calculating averages together with the pf_atom_*_merge() functions */
void pf_atom_summed_real_divide(t_pf_atom_summed *atom, real divisor) {
  t_pf_interaction_array_summed *ia;
  atom_id j;
  
  ia = &(atom->interactions);
  for (j = 0; j < ia->len; j++) {
    svdiv(divisor, ia->array[j].force);
  }
}
