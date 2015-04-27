/*
 * Functions related to scalar interactions between atoms.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types/simple.h"
#include "pf_array_scalar.h"
#include "gmx_fatal.h"
#include "gromacs/utility/smalloc.h"
#include "vec.h"

/* test for allocated to avoid creating a special function to just free arrays
 *
 * 2014-10-06 B.Doser @@ The array must not be allocated from scratch every md frame.
 * It's enough to zero a->len.
 */
static inline void pf_interaction_array_scalar_init(t_pf_interaction_array_scalar *a) {
  //if (a->allocated) sfree(a->array);
  //a->array = NULL;
  //a->allocated = 0;
  a->len = 0;
}

void pf_atom_scalar_init(t_pf_atom_scalar *atom) {
  pf_interaction_array_scalar_init(&atom->interactions);
}

t_pf_interaction_scalar *pf_lookup_interaction_scalar(t_pf_interaction_array_scalar *ia, int jjnr) {
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
void pf_interaction_array_scalar_append(t_pf_interaction_array_scalar *ia, atom_id jjnr, int type, real force) {
  //fprintf(stderr, "pf_interaction_array_scalar_append ia_summed=%p array=%p len=%d allocated=%d jjnr=%d force=%f\n", ia, ia->array, ia->len, ia->allocated, jjnr, force);
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
  ia->array[ia->len].type = type;
  ia->array[ia->len].force = force;
  ia->len++;
}

/* looks for the interaction in the array corresponding to type; if found,
 * the force is added to the existing value; if not, a new interaction is
 * appended to the array
 */
void pf_atom_scalar_add(t_pf_atom_scalar *atom, atom_id jjnr, int type, real force) {
  t_pf_interaction_scalar *p;

  //fprintf(stderr, "pf_atom_scalar_add: jjnr=%d, type=%d, force=%e\n", jjnr, type, force);
  p = pf_lookup_interaction_scalar(&atom->interactions, jjnr);
  if (p == NULL) {
    pf_interaction_array_scalar_append(&atom->interactions, jjnr, type, force);
  } else {
    p->force += force;
    p->type |= type;
  }
}

/* merge interactions of src into dst */
void pf_atom_scalar_merge(t_pf_atom_scalar *dst, t_pf_atom_scalar *src) {
  t_pf_interaction_array_scalar *ia;
  t_pf_interaction_scalar *i;
  atom_id j;

  if (src->nr != dst->nr)
    gmx_fatal(FARGS, "Mismatch of PF atom scalar number: %d vs %d\n", src->nr, dst->nr);
  ia = &(src->interactions);
  for (j = 0; j < ia->len; j++) {
    i = &(ia->array[j]);
    pf_atom_scalar_add(dst, i->jjnr, i->type, i->force);
  }
}

/* divides all forces by the given number; useful for calculating averages together with the pf_atom_*_merge() functions */
void pf_atom_scalar_real_divide(t_pf_atom_scalar *atom, real divisor) {
  t_pf_interaction_array_scalar *ia;
  atom_id j;

  ia = &(atom->interactions);
  for (j = 0; j < ia->len; j++) {
    ia->array[j].force /= divisor;
  }
}
