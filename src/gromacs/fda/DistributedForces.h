/*
 * DistributedForces.h
 *
 *  Created on: Oct 31, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_
#define SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_

#include <array>
#include <map>
#include <vector>
#include "ForceType.h"
#include "InteractionType.h"
#include "OnePair.h"
#include "ResultType.h"
#include "Vector.h"

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
  DistributedForces(ResultType result_type, OnePair one_pair, int syslen);

  /// Divide all scalar forces by the divisor
  void scalar_real_divide(real divisor);

  ///
  void summed_merge_to_scalar(const rvec *x, int Vector2Scalar);

private:

  friend class FDA;

  /// Result type
  ResultType result_type;

  /// Detailed or summed storage
  OnePair one_pair;

  /// Indexing table: real atom nr. to index in the pf array; this has length equal to the total nr. of atoms in system
  std::vector<int> sys2pf;

  /// Number of interaction nodes
  std::vector<int> syslen;

  /// Scalar values of forces
  std::map<int, std::map<int, real>> scalar;

  /// Vector values of forces
  std::map<int, std::map<int, Vector>> summed;

  /// Detailed values of forces
  std::map<int, std::map<int, std::array<Vector, InteractionType::num>>> detailed;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_ */
