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
#include "Vector2Scalar.h"

namespace fda {

struct Atom {};
struct Residue {};

template <class Base>
class FDABase
{
public:

  FDABase(ResultType result_type, OnePair one_pair, int syslen, std::string const& result_filename, bool no_end_zeros, Vector2Scalar v2s)
   : result_type(result_type),
	 distributed_forces(result_type, one_pair, syslen),
	 total_forces(),
	 result_file(result_filename),
	 no_end_zeros(no_end_zeros),
	 one_pair(one_pair),
	 v2s(v2s)
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

  void write_frame(rvec *x, int nsteps);

  void write_frame_detailed(rvec *x, bool bVector, int nsteps);

  void write_frame_summed(rvec *x, bool bVector, int nsteps);

  void write_frame_scalar(int nsteps);

  void sum_total_forces(rvec *x);

  void write_total_forces();

  void write_frame_atoms_scalar_compat(DistributedForces const& forces, FILE *f, int *framenr, gmx_bool ascii);

  /**
   * Main function for scalar time averages; saves data and decides when to write it out
   *
   * dealing with residues is more complicated because COMs have to be averaged over time;
   * it's not possible to average the atom positions and calculate COMs only once before
   * writing because for this fda->atoms would need to be initialized (to get the atom
   * number or to get the sys2ps mapping) which only happens when AtomBased is non-zero
   */
  void save_and_write_scalar_time_averages(rvec const *x, gmx_mtop_t *mtop);

  /**
   * Write scalar time averages; this is similar to pf_write_frame, except that time averages are used
   *
   * There are several cases:
   * steps == 0 - there was no frame saved at all or there are no frames saved since the last writing, so no writing needed
   * steps > 0, period == 0 - no frames written so far, but frames saved, so writing one frame
   * steps > 0, period != 0 - frames saved, so writing the last frame
   *
   * if this is called from pf_write_and_save... then steps is certainly > 0
   * period is certainly != 1, otherwise the averages functions are not called
   */
  void write_scalar_time_averages();

  /// The stress is the negative atom_vir value.
  void write_atom_virial_sum(FILE *f, tensor *atom_vir, int natoms);

  /// For von Mises no negative values are needed, since all items are squared.
  void write_atom_virial_sum_von_mises(FILE *f, tensor *atom_vir, int natoms);

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
