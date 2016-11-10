/*
 * DistributedForces.h
 *
 *  Created on: Oct 31, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_
#define SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_

#include <vector>
#include "DistributedForcesDetailed.h"
#include "DistributedForcesScalar.h"
#include "DistributedForcesVector.h"
#include "FDASettings.h"
#include "ForceType.h"
#include "ResultType.h"

namespace fda {

/**
 * Storage container for distributed forces
 * Same structure for atom and residue based distribution
 *
 * TODO: the bonded interaction type  might be possible to be reduced to a
 * single real vector (=eliminate jjnr) as the atoms interacting are known
 * at simulation start and do not change (also their order in the bond list
 * doesn't change). This would however require an extra step of looking up
 * the atom indices when writing out the force matrix. This would also require
 * a change in the way the atom structure is accessed, it makes no sense to
 * keep an array of t_pf_interaction items, but only an array of reals
 * representing forces.
 */
class DistributedForces
{
public:

  /// Constructor
  DistributedForces(ForceType force_type, FDASettings fda_settings);

  /// Divide all scalar forces by the divisor
  void scalar_real_divide(real divisor);

  ///
  void summed_merge_to_scalar(const rvec *x, int Vector2Scalar);

  ResultType get_result_type() const {
	if (force_type == ForceType::ATOMS) return fda_settings.atom_based_result_type;
	else if (force_type == ForceType::RESIDUES) return fda_settings.residue_based_result_type;
	else return ResultType::INVALID;
  }

  bool compatibility_mode() const {
    ResultType r = get_result_type();
    return r == ResultType::COMPAT_BIN or r == ResultType::COMPAT_ASCII;
  }

  bool stress_mode() const {
	ResultType r = get_result_type();
    return r == ResultType::PUNCTUAL_STRESS or r == ResultType::VIRIAL_STRESS or r == ResultType::VIRIAL_STRESS_VON_MISES;
  }

  bool PF_or_PS_mode() const {
	ResultType r = get_result_type();
	return r == ResultType::PAIRWISE_FORCES_VECTOR or r == ResultType::PAIRWISE_FORCES_SCALAR or r == ResultType::PUNCTUAL_STRESS;
  }

  bool VS_mode() const {
	ResultType r = get_result_type();
	return r == ResultType::VIRIAL_STRESS or r == ResultType::VIRIAL_STRESS_VON_MISES;
  }

private:

  /// Atom or residue based forces
  ForceType force_type;

  /// FDA settings
  FDASettings fda_settings;

  /// Indexing table: real atom nr. to index in the pf array; this has length equal to the total nr. of atoms in system
  std::vector<int> sys2pf;

  /// Scalar values of forces
  std::vector<DistributedForcesScalar> scalar;

  /// Vector values of forces
  std::vector<DistributedForcesVector> summed;

  /// Detailed values of forces
  std::vector<DistributedForcesDetailed> detailed;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_ */
