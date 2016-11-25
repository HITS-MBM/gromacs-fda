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
#include "gromacs/math/vectypes.h"
#include "PureInteractionType.h"
#include "OnePair.h"
#include "ResultType.h"
#include "Vector.h"
#include "Vector2Scalar.h"

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
  DistributedForces(ResultType result_type);

  /// Divide all scalar forces by the divisor
  void scalar_real_divide(real divisor);

  void summed_merge_to_scalar(const rvec *x, Vector2Scalar v2s);

  size_t size() const { return scalar.size(); }

private:

  friend class FDA;
  template <class Base> friend class FDABase;

  /// Result type
  ResultType result_type;

  /// Scalar values of forces
  std::map<int, std::map<int, real>> scalar;

  /// Vector values of forces
  std::map<int, std::map<int, Vector>> summed;

  /// Detailed values of forces
  //std::map<int, std::map<int, std::array<Vector, number_of_pure_interactions>>> detailed;
  std::map<int, std::map<int, std::array<Vector, static_cast<int>(PureInteractionType::NUMBER)>>> detailed;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_ */
