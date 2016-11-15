/*
 * FDASettings.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <set>
#include "FDASettings.h"

using namespace fda;

FDASettings::FDASettings(int nfile, const t_filenm fnm[], gmx_mtop_t *top_global)
 : atom_based_result_type(ResultType::NO),
   residue_based_result_type(ResultType::NO),
   one_pair(OnePair::DETAILED),
   vector_2_scalar(Vector2Scalar::NORM),
   residues_renumber(ResiduesRenumber::AUTO),
   no_end_zeros(false),
   syslen_atoms(top_global->natoms),
   syslen_residues(0),
   time_averaging_period(1),
   sys_in_group1(syslen_atoms, 0),
   sys_in_group2(syslen_atoms, 0),
   type(0)
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

  no_end_zeros = strcasecmp(get_estr(&ninp, &inp, "no_end_zeros", "no"), "no");

  if ((compatibility_mode(atom_based_result_type) or compatibility_mode(residue_based_result_type)) and Vector2Scalar != Vector2Scalar::NORM)
	gmx_fatal(FARGS, "When using compat mode, pf_vector2scalar should be set to norm.\n");

  std::stringstream(get_estr(&ninp, &inp, "type", "all")) >> type;
  std::cout << "Pairwise interactions selected: " << type << std::endl;

  if (type == InteractionType::NONE)
	gmx_fatal(FARGS, "No interactions selected, no sense to compute pairwise forces.\n");

  // Read group names
  std::string name_group1, name_group2;
  std::stringstream(get_estr(&ninp, &inp, "group1", "Protein")) >> name_group1;
  std::stringstream(get_estr(&ninp, &inp, "group2", "Protein")) >> name_group2;
  std::cout << "Pairwise forces for groups: " << name_group1 << " and " << name_group2 << std::endl;

  // Check if using compatibility mode, there should be only one group
  if (compatibility_mode(atom_based_result_type) or compatibility_mode(residue_based_result_type)) {
	if (name_group1 != name_group2)
	  gmx_fatal(FARGS, "When using compat mode, the two group names should the identical.\n");
	else
	  groupname = name_group1;
  }

  // Read group information
  char** group_names;
  t_blocka *groups = init_index(opt2fn("-pfn", nfile, fnm), &group_names);
  if(groups->nr == 0) gmx_fatal(FARGS, "No groups found in the indexfile.\n");

  // Set sys_in_group arrays
  for (int i = 0; i != groups->nr; ++i) {
    std::vector<int> group_atoms;
    if (group_names[i] == name_group1) {
      group_atoms = std::vector<int>(groups->a + groups->index[i], groups->a + groups->index[i + 1]);
      for (auto g : group_atoms) sys_in_group1[g] = 1;
    }
    if (group_names[i] == name_group2) {
      group_atoms = std::vector<int>(groups->a + groups->index[i], groups->a + groups->index[i + 1]);
      for (auto g : group_atoms) sys_in_group2[g] = 1;
    }
    if (PF_or_PS_mode(atom_based_result_type)) {
      for (auto g : group_atoms) sys2pf_atoms[g] = sys2pf_atoms.size();
    }
    if (PF_or_PS_mode(residue_based_result_type)) {
      std::vector<int> group_residues = groupatoms2residues(group_atoms);
      for (auto g : group_residues) sys2pf_residues[g] = sys2pf_residues.size();
    }
  }

  // Read time averaging period
  time_averaging_period = get_eint(&ninp, &inp, "time_averages_period", 1, wi);
  if (time_averaging_period < 0)
	gmx_fatal(FARGS, "Invalid value for time_averages_period: %d\n", time_averaging_period);

  // Check if groups are defined for PF/PF mode
  if (PF_or_PS_mode(atom_based_result_type) or PF_or_PS_mode(residue_based_result_type)) {
	if ((sys_in_group1 == NULL) || (sys_in_group2 == NULL))
	  gmx_fatal(FARGS, "No atoms in one or both groups.\n");
  }

  // Check that there is an index file
  if (!opt2bSet("-pfn", nfile, fnm))
	gmx_fatal(FARGS, "No index file (-pfn) for pairwise forces.\n");

  // Get output file names
  if (PF_or_PS_mode(atom_based_result_type)) {
	pf_atoms_alloc(OnePair, atom_based_forces, fda_settings.syslen_atoms, "atoms");
	if (atom_based_result_type == ResultType::PUNCTUAL_STRESS) {
		pf_per_atom_real_init(&(per_atom_real), fda_settings.syslen_atoms, 0.0);
		ofn_atoms = gmx_strdup(opt2fn("-psa", nfile, fnm));
	} else {
		ofn_atoms = gmx_strdup(opt2fn("-pfa", nfile, fnm));
	}
	if (!ofn_atoms)
	  gmx_fatal(FARGS, "No file for writing out atom-based values.\n");
  }

  if (PF_or_PS_mode(residue_based_result_type)) {
	pf_atoms_alloc(OnePair, residue_based_forces, syslen_residues, "residues");
	if (residue_based_result_type == ResultType::PUNCTUAL_STRESS) {
		pf_per_atom_real_init(&(per_residue_real), syslen_residues, 0.0);
		ofn_residues = gmx_strdup(opt2fn("-psr", nfile, fnm));
	} else {
		ofn_residues = gmx_strdup(opt2fn("-pfr", nfile, fnm));
	}
	if (!ofn_residues)
	  gmx_fatal(FARGS, "No file for writing out residue-based values.\n");
  }

  if (VS_mode(atom_based_result_type)) {
	snew(atom_vir, fda_settings.syslen_atoms);
	if (atom_based_result_type == ResultType::VIRIAL_STRESS)
		ofn_atoms = gmx_strdup(opt2fn("-vsa", nfile, fnm));
	if (atom_based_result_type == ResultType::VIRIAL_STRESS_VON_MISES)
		ofn_atoms = gmx_strdup(opt2fn("-vma", nfile, fnm));
	if (!ofn_atoms) gmx_fatal(FARGS, "No file for writing out virial stress.\n");
  }

  if ((stress_mode(atom_based_result_type) or stress_mode(residue_based_result_type)) and (one_pair != OnePair::SUMMED))
	gmx_fatal(FARGS, "Per atom data can only be computed from summed interactions.\n");

  fill_atom2residue(top_global);
}

std::vector<int> FDASettings::groupatoms2residues(std::vector<int> const& group_atoms) const
{
  std::set<int> group_residues;
  for (auto atom : group_atoms) group_residues.insert(atom_2_residue[atom]);
  return std::vector<int>(group_residues.begin(), group_residues.end());
}

void FDASettings::fill_atom2residue(gmx_mtop_t *top_global)
{
  int moltype_index, mol_index, atom_index, atom_global_index, resnrmax, i;
  t_atoms *atoms;
  t_atom *atom_info;
  gmx_molblock_t *mb;
  int *a2r_resnr, *a2r_renum;       /* atom 2 residue correspondence tables, both are filled, one will be used in the end */
  int *resnr2renum;                 /* residue nr. to renumbered nr. corespondence */
  int resnr;                            /* residue nr.; set to atoms->resinfo[].nr => type int */
  int renum;                        /* renumbered residue nr.; increased monotonically, so could theoretically be as large as the nr. of atoms => type int */
  gmx_bool bResnrCollision = FALSE;     /* TRUE if residue number collision */

  /* fill in the map between atom nr. and residue nr.; adapted from the code fragment in Data_Structures page on GROMACS website */
  snew(a2r_resnr, top_global->natoms);
  snew(a2r_renum, top_global->natoms);
  /* get the maximum number of residues in any molecule block;
   * any residue nr. obtained via atoms->resinfo[].nr is guaranteed to be smaller than this;
   * this will be used for residue nr. collision detection - initialized to -1,
   * will be set to the renumbered value on the first encounter, then a new renumbered value means a collision */
  resnrmax = gmx_mtop_maxresnr(top_global, 0);
  snew(resnr2renum, resnrmax + 1);
  //fprintf(stderr, "resnrmax=%d\n",resnrmax);
  for (i = 0; i <= resnrmax; i++)
    resnr2renum[i] = -1;
  atom_global_index = 0;
  renum = 0;
  for (moltype_index = 0; moltype_index < top_global->nmolblock; moltype_index++) {	//enum all molecule types
    mb = &top_global->molblock[moltype_index];
    for (mol_index = 0; mol_index < mb->nmol; mol_index++) {				//enum all instances of each molecule type
      atoms = &top_global->moltype[mb->type].atoms;
      for(atom_index = 0; atom_index < atoms->nr; atom_index++) {			//enum atoms
        atom_info=&atoms->atom[atom_index];
        resnr = atoms->resinfo[atom_info->resind].nr;
        renum = pf_get_global_residue_number(top_global, atom_global_index);
        //fprintf(stderr, "atom=%d, resnr=%d, renum=%d\n", atom_global_index, resnr, renum);
        //fprintf(stderr, "pair %d:%d\n", resnr, resnr2renum[resnr]);
        if ((resnr2renum[resnr] != renum) && (resnr2renum[resnr] != -1)) {
          bResnrCollision = TRUE;
          //fprintf(stderr, "Renumbering because %d:%d should be %d:%d\n", resnr, renum, resnr, resnr2renum[resnr]);
        }
        resnr2renum[resnr] = renum;
        a2r_resnr[atom_global_index] = resnr;
        a2r_renum[atom_global_index] = renum;
        atom_global_index++;
      }
    }
  } /* for (moltype_index = 0... */
  sfree(resnr2renum);

  /* renum is set to the residue number of the last atom, so should be largest value so far */
  switch (fda->ResiduesRenumber) {
    case PF_RESIDUE_RENUMBER_AUTO:
      if (bResnrCollision) {
        fprintf(stderr, "Residue number collision detected, residues will be renumbered.\n");
        fda->atom2residue = a2r_renum;
        fda->syslen_residues = renum + 1;
        sfree(a2r_resnr);
      }
      else {
        fprintf(stderr, "Residues will not be renumbered.\n");
        fda->atom2residue = a2r_resnr;
        fda->syslen_residues = resnrmax + 1;
        sfree(a2r_renum);
      }
      break;
    case PF_RESIDUE_RENUMBER_DO:
      fprintf(stderr, "Residues will be renumbered.\n");
      fda->atom2residue = a2r_renum;
      fda->syslen_residues = renum + 1;
      sfree(a2r_resnr);
      break;
    case PF_RESIDUE_RENUMBER_DONT:
      if (bResnrCollision)
        fprintf(stderr, "Residue number collision detected, residues will NOT be renumbered.\n");
      else
        fprintf(stderr, "Residues will not be renumbered.\n");
      fda->atom2residue = a2r_resnr;
      fda->syslen_residues = resnrmax + 1;
      sfree(a2r_renum);
      break;
  }
}
