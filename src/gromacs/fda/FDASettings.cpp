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
}


