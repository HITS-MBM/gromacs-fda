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
#include "OnePair.h"
#include "ResidueRenumber.h"
#include "ResultType.h"
#include "Vector2Scalar.h"

namespace fda {

/// Settings
struct FDASettings
{

  /// Default constructor
  FDASettings()
   : atom_based_result_type(ResultType::NO),
     residue_based_result_type(ResultType::NO),
	 one_pair(OnePair::DETAILED),
	 vector_2_scalar(Vector2Scalar::NORM),
	 residues_renumber(ResiduesRenumber::AUTO),
     syslen_atoms(0)
  {}

  /// Construction by input file
  FDASettings(int nfile, const t_filenm fnm[], gmx_mtop_t *top_global);

  /// Read atom/residue groups
  void read_group(const char *ndxfile, char *groupname, char **sys_in_g);

  bool compatibility_mode(ResultType const& r) const {
    return r == ResultType::COMPAT_BIN or r == ResultType::COMPAT_ASCII;
  }

  bool stress_mode(ResultType const& r) const {
    return r == ResultType::PUNCTUAL_STRESS or r == ResultType::VIRIAL_STRESS or r == ResultType::VIRIAL_STRESS_VON_MISES;
  }

  bool PF_or_PS_mode(ResultType const& r) const {
	return r == ResultType::PAIRWISE_FORCES_VECTOR or r == ResultType::PAIRWISE_FORCES_SCALAR or r == ResultType::PUNCTUAL_STRESS;
  }

  bool VS_mode(ResultType const& r) const {
	return r == ResultType::VIRIAL_STRESS or r == ResultType::VIRIAL_STRESS_VON_MISES;
  }

  /// ResultType for atom based forces
  ResultType atom_based_result_type;

  /// ResultType for residue based forces
  ResultType residue_based_result_type;

  /// OnePair defines the way the interactions between the same pair of atoms are stored
  OnePair one_pair;

  /// Define conversion from vector to scalar
  Vector2Scalar vector_2_scalar;

  /// detect/force residue renumbering
  ResiduesRenumber residues_renumber;

  /// Total number of atoms in the system.
  /// This is a local copy to avoid passing too many variables down the function call stack
  int syslen_atoms;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_FDASETTINGS_H_ */
