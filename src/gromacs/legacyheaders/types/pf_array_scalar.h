/*
 * Structures for storing scalar interactions between atoms.
 *
 * Copyright Bogdan Costescu 2011-2012
 */

#ifndef pf_types_array_scalar_h
#define pf_types_array_scalar_h

#include "gromacs/utility/real.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

/* type for a single interaction */
typedef struct {
  int jjnr;         /* atom index */
  real force;
  int type;             /* interaction type, will be an "or" of all interactions which have been added to force */
} t_pf_interaction_scalar;

/* type to describe an array of interactions */
typedef struct {
  t_pf_interaction_scalar *array;       /* array of interactions */
  int len;                              /* current number of items in array */
  int allocated;                        /* number of items allocated for array; always larger or equal to len */
} t_pf_interaction_array_scalar;

/* type for an atom; contains lists of all atoms this atom interacts with */
typedef struct {
  t_pf_interaction_array_scalar interactions;
  int nr;                           /* real atom nr. */
} t_pf_atom_scalar;

#endif  /* pf_types_array_scalar_h */
