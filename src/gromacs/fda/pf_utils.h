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

t_pf_int_list *pf_int_list_alloc(int len);
void pf_int_list_free(t_pf_int_list *p);
void pf_atoms_scalar_real_divide(t_pf_atoms *atoms, real divisor);
void pf_atom_summed_merge_to_scalar(t_pf_atom_summed *src, t_pf_atom_scalar *dst, const rvec *x, int Vector2Scalar);
void pf_atoms_summed_merge_to_scalar(t_pf_atoms *atoms, const rvec *x, int Vector2Scalar);
real pf_vector2signedscalar(const rvec v, const rvec xi, const rvec xj, int Vector2Scalar);
void pf_fill_sys2pf(int *sys2pf, int *len, t_pf_int_list *p);
char *pf_make_sys_in_group(int syslen, t_pf_int_list *p);
void pf_fill_atom2residue(FDA *fda, gmx_mtop_t *top_global);
t_pf_int_list *pf_group2atoms(int len, int *list);
t_pf_int_list *pf_groupatoms2residues(t_pf_int_list *atoms, FDA *fda);
void pf_check_sys_in_g(FDA *fda);
int pf_interactions_type_str2val(char *typestr);
char *pf_interactions_type_val2str(int type);
void pf_atoms_alloc(int OnePair, t_pf_atoms *atoms, int syslen, char *name);
void pf_atoms_init(int OnePair, t_pf_atoms *atoms);
void pf_atoms_scalar_alloc(t_pf_atoms *atoms, int syslen, char *name);
void pf_atoms_scalar_init(t_pf_atoms *atoms);
void pf_atoms_and_residues_init(FDA *fda);

#endif  /* pf_utils_h */
