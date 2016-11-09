/*
 * FDASettings.h
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_FDASETTINGS_H_
#define SRC_GROMACS_FDA_FDASETTINGS_H_

#include "gromacs/commandline/filenm.h"
#include "gromacs/topology/topology.h"
#include "ResultType.h"

namespace fda {

/// Settings
struct FDASettings
{

  /// Default constructor
  FDASettings()
   : atom_based_result_type(ResultType::NO),
     residue_based_result_type(ResultType::NO),
     syslen_atoms(0)
  {}

  /// Construction by input file
  FDASettings(int nfile, const t_filenm fnm[], gmx_mtop_t *top_global);

  /// ResultType for atom based forces
  ResultType atom_based_result_type;

  /// ResultType for residue based forces
  ResultType residue_based_result_type;

  /// Total number of atoms in the system.
  /// This is a local copy to avoid passing too many variables down the function call stack
  int syslen_atoms;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_FDASETTINGS_H_ */
