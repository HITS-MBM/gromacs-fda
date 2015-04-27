/*
 * Structures for storing detailed interactions between atoms.
 * 
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_types_array_detailed_h
#define pf_types_array_detailed_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "simple.h"

/* type for a single interaction */
typedef struct {
  atom_id jjnr; 		/* atom index */
  rvec force;
} t_pf_interaction_detailed;

/* type to describe an array of interactions */
typedef struct {
  t_pf_interaction_detailed *array;	/* array of interactions */
  int len;				/* current number of items in array */
  int allocated;			/* number of items allocated for array; always larger or equal to len */
} t_pf_interaction_array_detailed;

/* type to describe all interaction types for detailed mode */
typedef struct {
  t_pf_interaction_array_detailed coulomb;
  t_pf_interaction_array_detailed lj;
  t_pf_interaction_array_detailed nb14;
  t_pf_interaction_array_detailed bonds;
  t_pf_interaction_array_detailed angles;
  t_pf_interaction_array_detailed dihedrals;
  t_pf_interaction_array_detailed polar;
} t_pf_interactions_detailed;

/* type for an atom; contains lists of all atoms this atom interacts with */
typedef struct {
  t_pf_interactions_detailed interactions;
  atom_id nr;				/* real atom nr. */
} t_pf_atom_detailed;

#endif  /* pf_types_array_detailed_h */
