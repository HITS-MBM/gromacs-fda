/*
 * Functions related to summed interactions between atoms.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_array_summed_h
#define pf_array_summed_h

#include "types/pf_array_summed.h"

static inline void pf_interaction_array_summed_init(t_pf_interaction_array_summed *a);
void pf_atom_summed_init(t_pf_atom_summed *atom);
t_pf_interaction_summed *pf_lookup_interaction_summed(t_pf_interaction_array_summed *ia, int jjnr);
void pf_interaction_array_summed_append(t_pf_interaction_array_summed *ia, atom_id jjnr, int type, rvec force);
void pf_atom_summed_add(t_pf_atom_summed *atom, atom_id jjnr, int type, rvec force);
void pf_atom_summed_merge(t_pf_atom_summed *dst, t_pf_atom_summed *src);
void pf_atom_summed_real_divide(t_pf_atom_summed *atom, real divisor);

#endif  /* pf_array_summed_h */
