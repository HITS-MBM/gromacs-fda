/*
 * Structures for storing per-atom data.
 * 
 * Copyright Bogdan Costescu 2011-2012
 */

#ifndef pf_types_per_atom_h
#define pf_types_per_atom_h

#include "gromacs/utility/real.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

/* type for scalar interactions per atom that require only a real;
 * could be used for sum, min/max, etc.
 */
typedef struct {
  int len;          /* length of the below array, should be equal to the syslen */
  real *force;
} t_pf_per_atom_real;

/* type for scalar interactions per atom that require a real and an integer;
 * could be used for averages per atom:
 * the forces will be added up in force, then divided by interactions;
 * could also add something like standard deviation, but not sure what this could be shown graphically...
 */
typedef struct {
  int len;          /* length of the below arrays, should be equal to the syslen */
  real *force;
  int *interactions;
} t_pf_per_atom_real_int;

#endif  /* pf_types_per_atom_h */
