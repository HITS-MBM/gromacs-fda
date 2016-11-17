/*
 * Functions related to scalar interactions between atoms.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_array_scalar_h
#define pf_array_scalar_h

static inline void pf_interaction_array_scalar_init(t_pf_interaction_array_scalar *a);
void pf_atom_scalar_init(t_pf_atom_scalar *atom);
t_pf_interaction_scalar *pf_lookup_interaction_scalar(t_pf_interaction_array_scalar *ia, int jjnr);
void pf_interaction_array_scalar_append(t_pf_interaction_array_scalar *ia, int jjnr, int type, real force);
void pf_atom_scalar_add(t_pf_atom_scalar *atom, int jjnr, int type, real force);
void pf_atom_scalar_merge(t_pf_atom_scalar *dst, t_pf_atom_scalar *src);
void pf_atom_scalar_real_divide(t_pf_atom_scalar *atom, real divisor);

#endif  /* pf_array_scalar_h */
