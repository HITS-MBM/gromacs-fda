/*
 * Top level functions to add interactions.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_array_h
#define pf_array_h

#include "pf_array_detailed.h"
#include "pf_array_scalar.h"
#include "pf_array_summed.h"
#include "types/pf_array.h"

#ifdef __cplusplus
extern "C" {
#endif

void pf_atom_add_bonded(t_pf_global *pf_global, int i, int j, int type, rvec force);
void pf_atom_add_nonbonded_single(t_pf_global *pf_global, int i, int j, int type, real force, real dx, real dy, real dz);
void pf_atom_add_nonbonded(t_pf_global *pf_global, int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz);

void pf_atom_virial_add(t_pf_global *pf_global, int ai, tensor v, real s);
void pf_atom_virial_bond(t_pf_global *pf_global, int ai, int aj, real fbond, real dx, real dy, real dz);
void pf_atom_virial_angle(t_pf_global *pf_global, int ai, int aj, int ak, rvec r_ij, rvec r_kj, rvec f_i, rvec f_k);
void pf_atom_virial_dihedral(t_pf_global *pf_global, int i, int j, int k, int l, rvec f_i, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl);

#ifdef __cplusplus
}
#endif

#endif  /* pf_array_h */
