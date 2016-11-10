/*
 * FDASettings.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "FDASettings.h"

using namespace fda;

FDASettings::FDASettings(int nfile, const t_filenm fnm[], gmx_mtop_t *top_global)
 : atom_based_result_type(ResultType::NO),
   residue_based_result_type(ResultType::NO),
   one_pair(OnePair::DETAILED),
   vector_2_scalar(Vector2Scalar::NORM),
   residues_renumber(ResiduesRenumber::AUTO),
   syslen_atoms(top_global->natoms)
{
  // check for the pf configuration file (specified with -pfi option);
  // if it doesn't exist, return NULL to specify that no pf handling is done;
  // otherwise, check also for specification of the index file (-pfn)
  if (opt2bSet("-pfi", nfile, fnm)) {
	if (!opt2bSet("-pfn", nfile, fnm))
	  gmx_fatal(FARGS, "-pfi option (pairwise forces configuration) specified, an index file (-pfn) is also needed.\n");
  } else {
	gmx_fatal(FARGS, "No pairwise forces input file, no sense to compute pairwise forces.\n");
	return;
  }

  // Use GROMACS function to read lines of form key = value from input file
  const char *pf_file_in = opt2fn("-pfi", nfile, fnm);
  warninp_t wi = init_warning(FALSE, 0);
  int ninp;
  t_inpfile *inp = read_inpfile(pf_file_in, &ninp, wi);
  std::stringstream(get_estr(&ninp, &inp, "atombased", "no")) >> atom_based_result_type;
  std::stringstream(get_estr(&ninp, &inp, "residuebased", "no")) >> residue_based_result_type;

  // OnePair has to be initialized before the atoms/residues are initialized
  // because the data structures used for storing atoms/residues depend on it
  std::stringstream(get_estr(&ninp, &inp, "onepair", "detailed")) >> one_pair;

  std::stringstream(get_estr(&ninp, &inp, "residuesrenumber", "auto")) >> residues_renumber;
  std::cout << "ResidueRenumber: " << residues_renumber << std::endl;

  std::stringstream(get_estr(&ninp, &inp, "vector2scalar", "norm")) >> vector_2_scalar;
  std::cout << "Vector2Scalar: " << vector_2_scalar << std::endl;

  if ((compatibility_mode(atom_based_result_type) or compatibility_mode(residue_based_result_type)) and Vector2Scalar != Vector2Scalar::NORM)
	gmx_fatal(FARGS, "When using compat mode, pf_vector2scalar should be set to norm.\n");

  std::stringstream(get_estr(&ninp, &inp, "type", "all")) >> type;
  std::cout << "Pairwise interactions selected: " << type << std::endl;

  if (type == InteractionType::NONE)
	gmx_fatal(FARGS, "No interactions selected, no sense to compute pairwise forces.\n");

  // Obtain in g1name and g2name the names of the 2 groups
  std::string name_group1, name_group2;
  std::stringstream(get_estr(&ninp, &inp, "group1", "Protein")) >> name_group1;
  std::stringstream(get_estr(&ninp, &inp, "group2", "Protein")) >> name_group2;

  this->read_group(opt2fn("-pfn", nfile, fnm), name_group1, &sys_in_group1);
  this->read_group(opt2fn("-pfn", nfile, fnm), name_group2, &sys_in_group2);
}

void FDASettings::read_group(const char *ndxfile, char *groupname, char **sys_in_g)
{
  t_blocka *groups;
  char** groupnames;
  int i;
  t_pf_int_list *groupatoms, *groupresidues;
  int *residues_in_g;

  groups = init_index(ndxfile, &groupnames);
  if(groups->nr == 0)
    gmx_fatal(FARGS, "No groups found in the indexfile.\n");

  snew(residues_in_g, syslen_residues);

  for (i = 0; i < groups->nr; i++) {
    if (strcmp(groupnames[i], groupname) == 0) {
      //fprintf(stderr, "found group: %s\n", groupname);
      groupatoms = pf_group2atoms(groups->index[i + 1] - groups->index[i], &groups->a[groups->index[i]]);
      *sys_in_g = pf_make_sys_in_group(syslen_atoms, groupatoms);

      if (PF_or_PS_mode(AtomBased)) {
        pf_fill_sys2pf(atom_based_forces->sys2pf, &atom_based_forces->len, groupatoms);
//	for (k = 0; k < syslen_atoms; k++)
//	  fprintf(stderr, "sys2pf[%d]=%d\n", k, atoms->sys2pf[k]);
      }
      if (PF_or_PS_mode(ResidueBased)) {
        groupresidues = pf_groupatoms2residues(groupatoms, this);
        pf_fill_sys2pf(residue_based_forces->sys2pf, &residue_based_forces->len, groupresidues);
//	for (k = 0; k < syslen_residues; k++)
//	  fprintf(stderr, "sys2pf[%d]=%d\n", k, residues->sys2pf[k]);
        pf_int_list_free(groupresidues);
      }

      pf_int_list_free(groupatoms);
      break;	/* once the group is found, no reason to check against further group names */
    }
  }
}
