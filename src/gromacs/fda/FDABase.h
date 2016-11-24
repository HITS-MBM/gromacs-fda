/*
 * FDABase.h
 *
 *  Created on: Nov 22, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_FDABASE_H_
#define SRC_GROMACS_FDA_FDABASE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include "DistributedForces.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "OnePair.h"
#include "ResultType.h"
#include "Tensor.h"
#include "Vector2Scalar.h"

namespace fda {

/// Type for atom-based forces
/// Virial stress is only needed for atom-based forces
struct Atom
{
  Atom(bool VS_mode, int syslen)
   : virial_stress(VS_mode ? syslen : 0)
  {}

  /// Virial stress
  std::vector<Tensor> virial_stress;
};

/// Type for residue-based forces
struct Residue
{
  Residue(bool, int) {}
};

/**
 * @brief Managing atom- and residue-based operations
 *
 * Atom- and residue-based operations are almost identical except the dimension and indexing.
 * @FDA holds two instances: one for the atom-based and one for the residue-based forces.
 */
template <class Base>
class FDABase : Base
{
public:

  FDABase(ResultType result_type, OnePair one_pair, int syslen, std::string const& result_filename,
    bool no_end_zeros, Vector2Scalar v2s);

  bool compatibility_mode() const {
    return result_type == ResultType::COMPAT_BIN or
    	   result_type == ResultType::COMPAT_ASCII;
  }

  bool stress_mode() const {
    return result_type == ResultType::PUNCTUAL_STRESS or
    	   result_type == ResultType::VIRIAL_STRESS or
		   result_type == ResultType::VIRIAL_STRESS_VON_MISES;
  }

  bool PF_or_PS_mode() const {
	return result_type == ResultType::PAIRWISE_FORCES_VECTOR or
		   result_type == ResultType::PAIRWISE_FORCES_SCALAR or
		   result_type == ResultType::PUNCTUAL_STRESS;
  }

  bool VS_mode() const {
	return result_type == ResultType::VIRIAL_STRESS or
		   result_type == ResultType::VIRIAL_STRESS_VON_MISES;
  }

  void write_frame(rvec *x, int nsteps);

  void write_frame_detailed(rvec *x, bool bVector, int nsteps);

  void write_frame_summed(rvec *x, bool bVector, int nsteps);

  void write_frame_scalar(int nsteps);

  void sum_total_forces(rvec *x);

  void write_total_forces();

  void write_frame_atoms_compat();

  void write_frame_atoms_scalar_compat();

  void write_frame_atoms_summed_compat(rvec *x);

  /**
   * Writes a header as in original PF implementation;
   * as the original PF implementation calculated everything then wrote out everything,
   * nsteps was known at the time when the file was written; however, to make it work
   * when data is written as it comes (frame by frame), the header is written once when
   * the file is open, but with nsteps=1 and space for up to 8 digits; when the file is
   * closed, the header is written again, this time with the correct nsteps
   */
  void write_compat_header(int nsteps);

  /// The stress is the negative atom_vir value.
  void write_atom_virial_sum();

  /// For von Mises no negative values are needed, since all items are squared.
  void write_atom_virial_sum_von_mises();

private:

  friend class FDA;

  /// Result type
  ResultType result_type;

  /// Distributed forces
  DistributedForces distributed_forces;

  /// Total force per atom
  std::vector<real> total_forces;

  /// Result file
  std::ofstream result_file;

  /// If True, trim the line such that the zeros at the end are not written.
  /// if False (default), all per atom/residue data is written.
  bool no_end_zeros;

  /// OnePair defines the way the interactions between the same pair of atoms are stored
  OnePair one_pair;

  /// Define conversion from vector to scalar
  Vector2Scalar v2s;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_FDABASE_H_ */
