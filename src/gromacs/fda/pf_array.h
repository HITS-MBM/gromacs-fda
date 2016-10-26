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
#include "fda.h"

#ifdef __cplusplus
extern "C" {
#endif

void pf_atom_add_bonded_nocheck(t_pf_global *pf_global, int i, int j, int type, rvec force);

/**
 *  Add a particular type of nonbonded interaction for the kernels where only one type of interaction is computed;
 *  force is passed as scalar along with the distance vector (as dx, dy, dz) from which the vector force is
 *  computed, the same way it's done in the nonbonded kernels
 */
void pf_atom_add_nonbonded_single(t_pf_global *pf_global, int i, int j, int type, real force, real dx, real dy, real dz);

/** Add a nonbonded interaction for kernels where both Coulomb and LJ are computed;
 *  this is more efficient than calling the previous one twice because some of the tests are made only once;
 *  forces are passed as scalars along with the distance vector (as dx, dy, dz) from which the vector forces are
 *  computed, the same way it's done in the nonbonded kernels
 */
void pf_atom_add_nonbonded(t_pf_global *pf_global, int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz);

/**
 *  The atom virial can be expressed as a 6-real tensor, as it's symmetric.
 *  To avoid defining a new tensor type, the 9-real tensor is used instead.
 */
void pf_atom_virial_add(t_pf_global *pf_global, int ai, tensor v, real s);

/**
 *  Origin on j, but for 2 atoms it doesn't matter.
 */
void pf_atom_virial_bond(t_pf_global *pf_global, int ai, int aj, real fbond, real dx, real dy, real dz);

#ifdef __cplusplus
}
#endif

void pf_atom_add_bonded(t_pf_global *pf_global, int i, int j, int type, rvec force);

/**
 *  Translate to origin on the middle (j) atom:
 *  vir = ri*Fi + rj*Fj + rk*Fk
 *      = (ri-rj)*Fi + (rk-rj)*Fk
 *      = r_ij[dim1]*f_i[dim2] + r_kj[dim1]*f_k[dim2]
 */
void pf_atom_virial_angle(t_pf_global *pf_global, int ai, int aj, int ak, rvec r_ij, rvec r_kj, rvec f_i, rvec f_k);

/**
 *  Translate to origin on the second (j) atom:
 *  vir = ri*Fi + rj*Fj + rk*Fk + rl*Fl
 *      = (ri-rj)*Fi + (rk-rj)*Fk + (rl-rj)*Fl
 *      = (ri-rj)*Fi + (rk-rj)*Fk + ((rl-rk) + (rk-rj))*Fl
 *      = r_ij[dim1]*f_i[dim2] + r_kj[dim1]*f_k[dim2] + (r_kj-r_kl)[dim1]*f_l[dim2]
 */
void pf_atom_virial_dihedral(t_pf_global *pf_global, int i, int j, int k, int l, rvec f_i, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl);

#endif  /* pf_array_h */
