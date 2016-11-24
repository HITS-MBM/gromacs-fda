/*
 * Top level functions to add interactions.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_array_h
#define pf_array_h

#include "fda.h"

#ifdef __cplusplus
extern "C" {
#endif

void pf_atom_add_bonded_nocheck(struct FDA *fda, int i, int j, int type, rvec force);

/**
 *  Add a particular type of nonbonded interaction for the kernels where only one type of interaction is computed;
 *  force is passed as scalar along with the distance vector (as dx, dy, dz) from which the vector force is
 *  computed, the same way it's done in the nonbonded kernels
 */
void pf_atom_add_nonbonded_single(struct FDA *fda, int i, int j, int type, real force, real dx, real dy, real dz);

/** Add a nonbonded interaction for kernels where both Coulomb and LJ are computed;
 *  this is more efficient than calling the previous one twice because some of the tests are made only once;
 *  forces are passed as scalars along with the distance vector (as dx, dy, dz) from which the vector forces are
 *  computed, the same way it's done in the nonbonded kernels
 */
void pf_atom_add_nonbonded(struct FDA *fda, int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz);

/**
 *  The atom virial can be expressed as a 6-real tensor, as it's symmetric.
 *  To avoid defining a new tensor type, the 9-real tensor is used instead.
 */
void pf_atom_virial_add(struct FDA *fda, int ai, tensor v, real s);

/**
 *  Origin on j, but for 2 atoms it doesn't matter.
 */
void pf_atom_virial_bond(struct FDA *fda, int ai, int aj, real fbond, real dx, real dy, real dz);

#ifdef __cplusplus
}
#endif

#endif  /* pf_array_h */
