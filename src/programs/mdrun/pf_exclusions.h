/*
 * pf_exclusions.h
 *
 *  Created on: Sep 29, 2014
 *      Author: Bernd Doser, HITS gGmbH
 */

#ifndef PF_EXCLUSIONS_H_
#define PF_EXCLUSIONS_H_

#include "gromacs/fileio/filenm.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/utility/cstringutil.h"

struct gmx_groups_t;
struct gmx_mtop_t;
struct t_blocka;

#define PF_GROUP_IDX_FDA1  0
#define PF_GROUP_IDX_FDA2  1
#define PF_GROUP_IDX_FDA12 2
#define PF_GROUP_IDX_REST  3
#define PF_GROUP_DIM       4
#define PF_GROUP_DIM2      PF_GROUP_DIM*PF_GROUP_DIM

static int FDA_egp_flags[PF_GROUP_DIM2] = {
    1, 0, 0, 1,
    0, 1, 0, 1,
    0, 0, 0, 1,
    1, 1, 1, 1
};

#define FDA_PRINT_DEBUG_OFF
#define FDA_BONDEXCL_PRINT_DEBUG_OFF

typedef struct {
	char          g1name[STRLEN]; /* groupname of FDA group 1 (defined in pfi-file) */
	char          g2name[STRLEN]; /* groupname of FDA group 2 (defined in pfi-file) */
    int           g1idx; /* FDA groups index 1 (defined in pfi-file) */
    int           g2idx; /* FDA groups index 2 (defined in pfi-file) */
	t_blocka*     groups; /* groups defined in pfn-file */
	char**        groupnames;/* groupnames defined in pfn-file */
	int           type; /* interaction type */
	char          FDA_nonbonded_exclusion_on; /* switch FDA nonbonded exclusions on/off (default: on) */
	char          FDA_bonded_exclusion_on; /* switch FDA bonded exclusions on/off (default: on) */
} t_pf_global_data;

//! Global data structure
//! Parts are redundant to t_pf_global, which is located in t_forcerec.
//! The function doing the nonbonded and bonded exclusions are called before the t_forcerec is initialized.
//! In the current implementation a second global object (real global instead of pf_global)
//! is used only for the data needed by the exclusion code.
//! For the future both global object should be merged.
t_pf_global_data pf_global_data;

#ifdef __cplusplus
extern "C" {
#endif

//! Initialization of global data
void pf_global_data_init(int nfile, const t_filenm fnm[]);

//! Append group to energy groups, returns the position index
int pf_add_name_to_energygrp(char* name, gmx_groups_t* groups);

//! Return position of group in energy groups array, if not found exit with fatal error
int pf_get_index_in_energygrp(char const* name, gmx_groups_t const* groups);

//! FDA groups must not be defined over complete charge groups.
//! This group redefine the energy group array with respect to the charge groups.
void pf_respect_charge_groups(unsigned char* array, gmx_mtop_t const* mtop);

//! Print exclusion table as matrix
void pf_print_exclusion_table(int* egp_flags, int dim);

//! Main routine for FDA exclusions
void pf_modify_energy_group_exclusions(gmx_mtop_t *mtop, t_inputrec *inputrec);

#ifdef __cplusplus
}
#endif

#endif /* PF_EXCLUSIONS_H_ */
