/*
 * Functions for per-atom interactions.
 *
 * Copyright Bogdan Costescu 2011-2013
 */

#ifndef pf_per_atom_h
#define pf_per_atom_h

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "types/pf_per_atom.h"
#include "types/pf_array_summed.h"

#ifdef __cplusplus
extern "C" {
#endif

void pf_per_atom_real_set(t_pf_per_atom_real *per_atom_real, real val);
void pf_per_atom_real_init(t_pf_per_atom_real **per_atom_real, int len, real val);
void pf_per_atom_real_int_set(t_pf_per_atom_real_int *per_atom_real_int, real val_real, int val_int);
void pf_per_atom_real_int_init(t_pf_per_atom_real_int **per_atom_real_int, int len, real val_real, int val_int);
void pf_per_atom_real_write_frame(FILE *f, real *force, int len, gmx_bool no_end_zeros);

void pf_per_atom_sum(t_pf_per_atom_real *per_atom_real, t_pf_atom_summed *atoms, int atoms_len, const rvec *x, int Vector2Scalar);
void pf_per_atom_average(t_pf_per_atom_real_int *per_atom_average, t_pf_atom_summed *atoms, int atoms_len, const rvec *x, int Vector2Scalar);
void pf_per_atom_minmax(t_pf_per_atom_real *per_atom_real, t_pf_atom_summed *atoms, int atoms_len, gmx_bool findmax, const rvec *x, int Vector2Scalar);

void pf_write_atom_virial_sum(FILE *f, tensor *atom_vir, int natoms);
void pf_write_atom_virial_sum_von_mises(FILE *f, tensor *atom_vir, int natoms);

#ifdef __cplusplus
}
#endif

#endif  /* pf_per_atom_h */
