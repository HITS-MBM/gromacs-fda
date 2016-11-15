/*
 * Functions which act on types/pf_array structures.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_utils_h
#define pf_utils_h

#include "fda.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/topology/topology.h"

real pf_vector2signedscalar(const rvec v, const rvec xi, const rvec xj, int Vector2Scalar);

int pf_interactions_type_str2val(char *typestr);
char *pf_interactions_type_val2str(int type);

#endif  /* pf_utils_h */
