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
#include "OnePair.h"
#include "ResultType.h"

namespace fda {

struct Atom {};
struct Residue {};

template <class Base>
class FDABase
{
public:

  FDABase(ResultType result_type, OnePair one_pair, int syslen, std::string const& result_filename, bool no_end_zeros)
   : result_type(result_type),
	 distributed_forces(result_type, one_pair, syslen),
	 total_forces(),
	 result_file(result_filename),
	 nsteps(0),
	 no_end_zeros(no_end_zeros)
  {
    if (PF_or_PS_mode()) {
	  if (result_type == ResultType::PUNCTUAL_STRESS) {
        total_forces.resize(syslen, 0.0);
      }
    }
  }

  bool compatibility_mode() const {
    return result_type == ResultType::COMPAT_BIN or result_type == ResultType::COMPAT_ASCII;
  }

  bool stress_mode() const {
    return result_type == ResultType::PUNCTUAL_STRESS or result_type == ResultType::VIRIAL_STRESS or result_type == ResultType::VIRIAL_STRESS_VON_MISES;
  }

  bool PF_or_PS_mode() const {
	return result_type == ResultType::PAIRWISE_FORCES_VECTOR or result_type == ResultType::PAIRWISE_FORCES_SCALAR or result_type == ResultType::PUNCTUAL_STRESS;
  }

  bool VS_mode() const {
	return result_type == ResultType::VIRIAL_STRESS or result_type == ResultType::VIRIAL_STRESS_VON_MISES;
  }

  void real_write_frame();

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

  /// Number of steps for output in compatibility mode
  /// Incremented for each step written during run, also used to write the total number of steps at the end
  int nsteps;

  /// If True, trim the line such that the zeros at the end are not written.
  /// if False (default), all per atom/residue data is written.
  bool no_end_zeros;

};

template <class Base>
void FDABase<Base>::real_write_frame()
{
  int j = total_forces.size();
  // Detect the last non-zero item
  if (no_end_zeros) {
    for (; j > 0; --j)
      if (total_forces[j - 1] != 0.0)
        break;
  }

  // j holds the index of first zero item or the length of force
  bool first_on_line = true;
  for (int i = 0; i < j; ++i) {
    if (first_on_line) {
      result_file << total_forces[i];
      first_on_line = false;
    } else {
      result_file << " " << total_forces[i];
    }
  }
  result_file << std::endl;
}

} // namespace fda

#endif /* SRC_GROMACS_FDA_FDABASE_H_ */
