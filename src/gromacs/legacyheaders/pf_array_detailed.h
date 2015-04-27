/*
 * Functions related to detailed interactions between atoms.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_array_detailed_h
#define pf_array_detailed_h

#include "types/pf_array_detailed.h"

static inline void pf_interaction_array_detailed_init(t_pf_interaction_array_detailed *a);
void pf_interactions_detailed_init(t_pf_interactions_detailed *interactions);
void pf_atom_detailed_init(t_pf_atom_detailed *atom);
t_pf_interaction_array_detailed *pf_interaction_array_by_type(t_pf_interactions_detailed *interactions, int type);
t_pf_interaction_detailed *pf_lookup_interaction_detailed(t_pf_interaction_array_detailed *ia, int jjnr);
void pf_interaction_array_detailed_append(t_pf_interaction_array_detailed *ia, atom_id jjnr, rvec force);
void pf_atom_detailed_add(t_pf_atom_detailed *atom, atom_id jjnr, int type, rvec force);
void pf_atom_detailed_merge(t_pf_atom_detailed *dst, t_pf_atom_detailed *src);
void pf_atom_detailed_real_divide(t_pf_atom_detailed *atom, real divisor);

#endif  /* pf_array_detailed_h */
