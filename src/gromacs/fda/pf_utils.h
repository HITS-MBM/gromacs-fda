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

void pf_fill_sys2pf(int *sys2pf, int *len, t_pf_int_list *p);
char *pf_make_sys_in_group(int syslen, t_pf_int_list *p);
t_pf_int_list *pf_group2atoms(int len, int *list);
t_pf_int_list *pf_groupatoms2residues(t_pf_int_list *atoms, FDA *fda);

void pf_check_sys_in_g(FDA *fda);

int pf_interactions_type_str2val(char *typestr);
char *pf_interactions_type_val2str(int type);

#endif  /* pf_utils_h */
