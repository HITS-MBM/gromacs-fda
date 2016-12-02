/*
 * pf_exclusions.c
 *
 *  Created on: Sep 29, 2014
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstring>
#include "InteractionType.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "pf_exclusions.h"

void fda_data_init(int nfile, const t_filenm fnm[])
{
	// Open pfi-file
	const char *pf_file_in = opt2fn("-pfi", nfile, fnm);
	warninp_t wi = init_warning(FALSE,0);
	int ninp;
	const char *tmp;
	t_inpfile *inp = read_inpfile(pf_file_in, &ninp, wi);

	// Get group names defined in pfi-file
	STYPE("group1", fda_data.g1name, "Protein");
	STYPE("group2", fda_data.g2name, "Protein");

	#ifdef FDA_PRINT_DEBUG_ON
		printf("=== DEBUG === g1name = %s\n", fda_data.g1name);
		printf("=== DEBUG === g2name = %s\n", fda_data.g2name);
	#endif

	// Get atom indices of groups defined in pfn-file
	fda_data.groups = init_index(opt2fn("-pfn", nfile, fnm), &fda_data.groupnames);

	// Get group number
	int i;
	for (i = 0; i < fda_data.groups->nr; ++i) {
		if (gmx_strcasecmp(fda_data.groupnames[i], fda_data.g1name) == 0) fda_data.g1idx = i;
		if (gmx_strcasecmp(fda_data.groupnames[i], fda_data.g2name) == 0) fda_data.g2idx = i;
	}

	#ifdef FDA_PRINT_DEBUG_ON
		printf("=== DEBUG === g1idx = %i\n", fda_data.g1idx);
		printf("=== DEBUG === g2idx = %i\n", fda_data.g2idx);
	#endif

	// Get interaction type
    char tmpstr[STRLEN];
    STYPE("type", tmpstr, "all");
    // TODO:
    //fda_data.type = pf_interactions_type_str2val(tmpstr);

	STYPE("energy_grp_exclusion", tmpstr, "yes");
	fda_data.FDA_nonbonded_exclusion_on = gmx_strcasecmp(tmpstr, "no") != 0;

	STYPE("bonded_exclusion", tmpstr, "yes");
	fda_data.FDA_bonded_exclusion_on = gmx_strcasecmp(tmpstr, "no") != 0;
}

int pf_add_name_to_energygrp(char* name, gmx_groups_t* groups)
{
	int index = groups->ngrpname;
	if (index == 255)
		gmx_fatal(FARGS, "FDA error: Limit of energy groups (256) exceeded.");
	groups->ngrpname += 1;
	srenew(groups->grpname, groups->ngrpname);
	char **tmp1;
	snew(tmp1,1);
	tmp1[0] = strdup(&name[0]);
	groups->grpname[index] = tmp1;
	return index;
}

int pf_get_index_in_energygrp(char const* name, gmx_groups_t const* groups)
{
	int i, index = -1;
	for (i = 0; i < groups->ngrpname; ++i) {
		if (gmx_strcasecmp(*groups->grpname[i], name) == 0) {
			index = i;
			break;
		}
	}

	if (index == -1)
		gmx_fatal(FARGS, "FDA error: group \"rest\" not found.");
	return index;
}

void pf_respect_charge_groups(unsigned char* array, gmx_mtop_t const* mtop)
{
	bool found_FDA1, found_FDA2, found_REST;
	int i, j, k, l, length;
	int nbMolecules = 0;
	int molTypeIndex = 0;
	int nbChargeGroups = 0;
	int *chargeGroupsIndex = NULL;
	unsigned char* p = array;

	for (i = 0; i < mtop->nmolblock; ++i) {

		nbMolecules = mtop->molblock[i].nmol;
		molTypeIndex = mtop->molblock[i].type;
		nbChargeGroups = mtop->moltype[molTypeIndex].cgs.nr;
		chargeGroupsIndex = mtop->moltype[molTypeIndex].cgs.index;

		for (j = 0; j < nbMolecules; ++j) {
			for (k = 0; k < nbChargeGroups; ++k) {
               length = chargeGroupsIndex[k+1] - chargeGroupsIndex[k];
               found_FDA1 = false; found_FDA2 = false; found_REST = false;
               for (l = 0; l < length; ++l) {
            	   if (p[l] == PF_GROUP_IDX_FDA1) found_FDA1 = true;
            	   if (p[l] == PF_GROUP_IDX_FDA2) found_FDA2 = true;
            	   if (p[l] == PF_GROUP_IDX_REST) found_REST = true;
               }
               if (found_FDA1 && found_FDA2) for (l = 0; l < length; ++l) p[l] = PF_GROUP_IDX_FDA12;
               else if (found_FDA1 && found_REST) for (l = 0; l < length; ++l) p[l] = PF_GROUP_IDX_FDA1;
               else if (found_FDA2 && found_REST) for (l = 0; l < length; ++l) p[l] = PF_GROUP_IDX_FDA2;
               p += length;
			}
		}
	}
}

void pf_print_exclusion_table(int* egp_flags, int dim)
{
	int i, j;
	for (i = 0; i < dim; ++i) {
		for (j = 0; j < dim; ++j) {
			printf("%i ", egp_flags[i*dim + j]);
		}
		printf("\n");
	}
}

void pf_modify_energy_group_exclusions(gmx_mtop_t *mtop, t_inputrec *inputrec)
{
	if (!fda_data.FDA_nonbonded_exclusion_on) {
		printf("WARNING: FDA energy group exclusion is set off.\n");
		return;
	}

	// If no nonbonded interactions are needed we simply exclude all energy groups
	if (!(fda_data.type & static_cast<int>(fda::InteractionType::NONBONDED))) {

		int i;
		for (i = 0; i < inputrec->opts.ngener; ++i) inputrec->opts.egp_flags[i] = 1;

	// If only one energy group is defined, e.g. the automatic generated group 'rest',
	// the FDA groups can easily added to the energy groups and energy group exclusions.
	// There must only be checked that no charge group will be split.
    } else if (inputrec->opts.ngener == 1) {

		#ifdef FDA_PRINT_DEBUG_ON
			printf("FDA: No energy groups are defined.\n");
		#endif

		mtop->groups.ngrpnr[egcENER] = mtop->natoms;
		snew(mtop->groups.grpnr[egcENER],mtop->natoms);

		// First set all atoms to energy group "rest" ...
		int i;
		for (i = 0; i < mtop->natoms; ++i) {
			mtop->groups.grpnr[egcENER][i] = PF_GROUP_IDX_REST;
		}

		// ... and then overwrite the defined atoms with FDA group 1
		for (i = fda_data.groups->index[fda_data.g1idx]; i < fda_data.groups->index[fda_data.g1idx+1]; ++i) {
			mtop->groups.grpnr[egcENER][fda_data.groups->a[i]] = PF_GROUP_IDX_FDA1;
		}

		// ... and FDA group 2
		for (i = fda_data.groups->index[fda_data.g2idx]; i < fda_data.groups->index[fda_data.g2idx+1]; ++i) {
			if (mtop->groups.grpnr[egcENER][fda_data.groups->a[i]] == PF_GROUP_IDX_FDA1)
				mtop->groups.grpnr[egcENER][fda_data.groups->a[i]] = PF_GROUP_IDX_FDA12;
			else
				mtop->groups.grpnr[egcENER][fda_data.groups->a[i]] = PF_GROUP_IDX_FDA2;
		}

		pf_respect_charge_groups(mtop->groups.grpnr[egcENER],mtop);

		// Search FDA group names in mtop->groups (tpr-file)
		int mtop_g1idx = pf_add_name_to_energygrp("FDA1", &mtop->groups);
		int mtop_g2idx = pf_add_name_to_energygrp("FDA2", &mtop->groups);
		int mtop_g3idx = pf_add_name_to_energygrp("FDA12", &mtop->groups);
		int mtop_rest_idx = pf_get_index_in_energygrp("rest", &mtop->groups);

		#ifdef FDA_PRINT_DEBUG_ON
			printf("=== DEBUG === mtop_g1idx = %i\n", mtop_g1idx);
			printf("=== DEBUG === mtop_g2idx = %i\n", mtop_g2idx);
			printf("=== DEBUG === mtop_g3idx = %i\n", mtop_g3idx);
			printf("=== DEBUG === mtop_rest_idx = %i\n", mtop_rest_idx);
		#endif

		// Lookup table energy group index to group index
		mtop->groups.grps[egcENER].nr = PF_GROUP_DIM;
		snew(mtop->groups.grps[egcENER].nm_ind, PF_GROUP_DIM);
		mtop->groups.grps[egcENER].nm_ind[0] = mtop_g1idx;
		mtop->groups.grps[egcENER].nm_ind[1] = mtop_g2idx;
		mtop->groups.grps[egcENER].nm_ind[1] = mtop_g3idx;
		mtop->groups.grps[egcENER].nm_ind[2] = mtop_rest_idx;

		// Write egp_flags table
		inputrec->opts.ngener = PF_GROUP_DIM;
		snew(inputrec->opts.egp_flags,PF_GROUP_DIM2);
		for (i = 0; i < PF_GROUP_DIM2; ++i) inputrec->opts.egp_flags[i] = FDA_egp_flags[i];

	} else {

		#ifdef FDA_PRINT_DEBUG_ON
			printf("FDA: %i energy groups are defined.\n", inputrec->opts.ngener);
		#endif

		unsigned char* FDA_eg;
		snew(FDA_eg,mtop->natoms);

		// First set all atoms to energy group "rest" ...
		int i;
		for (i = 0; i < mtop->natoms; ++i) {
			FDA_eg[i] = PF_GROUP_IDX_REST;
		}

		// ... and then overwrite the defined atoms with FDA group 1
		for (i = fda_data.groups->index[fda_data.g1idx]; i < fda_data.groups->index[fda_data.g1idx+1]; ++i) {
			FDA_eg[fda_data.groups->a[i]] = PF_GROUP_IDX_FDA1;
		}

		// ... and FDA group 2
		for (i = fda_data.groups->index[fda_data.g2idx]; i < fda_data.groups->index[fda_data.g2idx+1]; ++i) {
			FDA_eg[fda_data.groups->a[i]] = PF_GROUP_IDX_FDA2;
		}

		pf_respect_charge_groups(FDA_eg, mtop);

		// new values with prefix c_ for combined groups
		int c_ngener = PF_GROUP_DIM * inputrec->opts.ngener;
		unsigned char* c_eg;
		snew(c_eg, mtop->natoms);

		// Generate new energy group indices as combination of original energy groups and FDA groups
		for (i = 0; i < mtop->natoms; ++i) {
			c_eg[i] = mtop->groups.grpnr[egcENER][i] * PF_GROUP_DIM + FDA_eg[i];
		}

		int* c_egp_flags;
		snew(c_egp_flags,c_ngener * c_ngener);

		#ifdef FDA_PRINT_DEBUG_ON
			printf("=== DEBUG === inputrec->opts.egp_flags\n");
			pf_print_exclusion_table(inputrec->opts.egp_flags, inputrec->opts.ngener);
		#endif

		#ifdef FDA_PRINT_DEBUG_ON
			printf("=== DEBUG === FDA_egp_flags\n");
			pf_print_exclusion_table(FDA_egp_flags, PF_GROUP_DIM);
		#endif

		// Build new exclusion table
		int j,k,l;
		for (i = 0; i < inputrec->opts.ngener; ++i) {
			for (j = 0; j < inputrec->opts.ngener; ++j) {
				for (k = 0; k < PF_GROUP_DIM; ++k) {
					for (l = 0; l < PF_GROUP_DIM; ++l) {
						c_egp_flags[(i*PF_GROUP_DIM+k)*c_ngener + (j*PF_GROUP_DIM+l)] =
							inputrec->opts.egp_flags[i*inputrec->opts.ngener+j] || FDA_egp_flags[k*PF_GROUP_DIM+l];
					}
				}
			}
		}

		// Determine start index of energy groups
		int startIdx = UCHAR_MAX;
		for (i = 0; i < mtop->groups.grps[egcENER].nr; ++i) {
			startIdx = std::min(startIdx, mtop->groups.grps[egcENER].nm_ind[i]);
		}

		#ifdef FDA_PRINT_DEBUG_ON
			printf("=== DEBUG === startIdx = %i\n", startIdx);
		#endif

		// Rename existing groups
		char buffer[15];
		for (i = 0; i < inputrec->opts.ngener; ++i) {
			sprintf(buffer, "FDA%d", i);
			char **tmp;
			snew(tmp,1);
			tmp[0] = strdup(&buffer[0]);
			mtop->groups.grpname[i+startIdx] = tmp;
		}

		// Add additional group names
		for (i = 0; i < inputrec->opts.ngener * (PF_GROUP_DIM - 1); ++i) {
			sprintf(buffer, "FDA%d", i + inputrec->opts.ngener);
			pf_add_name_to_energygrp(buffer, &mtop->groups);
		}

		// Update lookup table energy group index to group index
		mtop->groups.grps[egcENER].nr = c_ngener;
		srenew(mtop->groups.grps[egcENER].nm_ind, c_ngener);
		for (i = 0; i < c_ngener; ++i) {
			mtop->groups.grps[egcENER].nm_ind[i] = i + startIdx;
		}

		// Set new energy group array
		sfree(mtop->groups.grpnr[egcENER]);
		mtop->groups.grpnr[egcENER] = c_eg;

		// Set new egp_flags table
		inputrec->opts.ngener = c_ngener;
		sfree(inputrec->opts.egp_flags);
		inputrec->opts.egp_flags = c_egp_flags;

		sfree(FDA_eg);
	}

#ifdef FDA_PRINT_DEBUG_ON

	printf("=== DEBUG === Print final arrays:\n");

	printf("=== DEBUG === mtop->groups.ngrpname = %i\n", mtop->groups.ngrpname);

	int i;
	for (i = 0; i < mtop->groups.ngrpname; ++i) {
		printf("=== DEBUG === mtop->groups.grpname[%i] = %s\n", i, *mtop->groups.grpname[i]);
	}

	printf("=== DEBUG === mtop->groups.ngrpnr[egcENER] = %i\n", mtop->groups.ngrpnr[egcENER]);

	for (i = 0; i < mtop->natoms; ++i) {
		printf("=== DEBUG === ggrpnr(&mtop->groups, egcENER, %i) = %i\n", i, (unsigned int)ggrpnr(&mtop->groups, egcENER, i));
	}

	printf("=== DEBUG === inputrec->opts->ngener = %i\n", inputrec->opts.ngener);
	pf_print_exclusion_table(inputrec->opts.egp_flags, inputrec->opts.ngener);

	printf("=== DEBUG === mtop->nmoltype = %i\n", mtop->nmoltype);

	int j;
	for (i = 0; i < mtop->nmoltype; ++i) {
		printf("=== DEBUG === mtop->moltype[%i].cgs.nr = %i\n", i, mtop->moltype[i].cgs.nr);

		for (j = 0; j < mtop->moltype[i].cgs.nr + 1; ++j) {
			printf("=== DEBUG === mtop->moltype[%i].cgs.index[%i] = %i\n",
				i, j, mtop->moltype[i].cgs.index[j]);
		}
	}

	printf("=== DEBUG === mtop->nmolblock = %i\n", mtop->nmolblock);

	for (i = 0; i < mtop->nmolblock; ++i) {
		printf("=== DEBUG === mtop->molblock[%i].type = %i\n", i, mtop->molblock[i].type);
		printf("=== DEBUG === mtop->molblock[%i].nmol = %i\n", i, mtop->molblock[i].nmol);
	}

#endif

}
