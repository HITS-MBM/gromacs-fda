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
	 no_end_zeros(false),
     syslen_atoms(0),
     syslen_residues(0),
	 time_averaging_period(1),
     sys_in_g1(nullptr),
     sys_in_g2(nullptr),
     type(0)
  {}

  /// Construction by input file
  FDASettings(int nfile, const t_filenm fnm[], gmx_mtop_t *top_global);

  /// Read atom/residue groups
  void read_group(const char *ndxfile, char *groupname, char **sys_in_g);

  /// Returns true if atoms i and j are in fda groups
  gmx_bool atoms_in_groups(int i, int j) const {
	return ((sys_in_g1[i] && sys_in_g2[j]) || (sys_in_g1[j] && sys_in_g2[i]));
  }

  void fill_atom2residue(gmx_mtop_t *top_global);

  int get_atom2residue(int i) const { return atom_2_residue[i]; }

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

  /// If True, trim the line such that the zeros at the end are not written.
  /// if False (default), all per atom/residue data is written.
  gmx_bool no_end_zeros;

  /// Total number of atoms in the system.
  /// This is a local copy to avoid passing too many variables down the function call stack
  int syslen_atoms;

  /// Maximum of residue nr. + 1; residue nr. doesn't have to be continuous, there can be gaps
  int syslen_residues;

  /// Number of steps to average before writing.
  /// If 1 (default), no averaging is done.
  /// If 0 averaging is done over all steps so only one frame is written at the end.
  int time_averaging_period;

  /// Output file name for atoms if AtomBased is non-zero
  std::string ofn_atoms;

  /// Output file name for residues if ResidueBased is non-zero
  std::string ofn_residues;

  /// If 0 if atom not in group1, if 1 if atom in group1, length of syslen_atoms; always allocated
  char *sys_in_g1;

  /// If 0 if atom not in group2, if 1 if atom in group2, length of syslen_atoms; always allocated
  char *sys_in_g2;

  /// Name of group for output in compatibility mode
  std::string groupname;

  /// Interaction types that are interesting, set based on input file; functions are supposed to test against this before calculating/storing data
  int type;

  /// Stores the residue number for each atom; array of length syslen; only initialized if ResidueBased is non-zero
  std::vector<int> atom_2_residue;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_FDASETTINGS_H_ */
