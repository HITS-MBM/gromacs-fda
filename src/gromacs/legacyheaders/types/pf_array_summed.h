/*
 * Structures for storing summed interactions between atoms.
 * 
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_types_array_summed_h
#define pf_types_array_summed_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "simple.h"

/* type for a single interaction */
typedef struct {
  atom_id jjnr;		/* atom index */
  rvec force;
  int type;		/* interaction type, will be an "or" of all interactions which have been added to force */
} t_pf_interaction_summed;

/* type to describe an array of interactions */
typedef struct {
  t_pf_interaction_summed *array;	/* array of interactions */
  int len;				/* current number of items in array */
  int allocated;			/* number of items allocated for array; always larger or equal to len */
} t_pf_interaction_array_summed;

/* type for an atom; contains lists of all atoms this atom interacts with */
typedef struct {
  t_pf_interaction_array_summed interactions;
  atom_id nr;				/* real atom nr. */
} t_pf_atom_summed;

#endif  /* pf_types_array_summed_h */
