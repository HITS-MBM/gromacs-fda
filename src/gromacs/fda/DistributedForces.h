/*
 * DistributedForces.h
 *
 *  Created on: Oct 31, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_
#define SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_

#include <stddef.h>
#include <map>
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"
#include "DetailedForce.h"
#include "Force.h"
#include "Vector.h"
#include "Vector2Scalar.h"

/// Forwarding needed for friend declaration
class FDA;

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

  /// Divide all scalar forces by the divisor
  void scalar_real_divide(real divisor);

  void summed_merge_to_scalar(const rvec *x, Vector2Scalar v2s);

  size_t size() const { return scalar.size(); }

private:

  friend class ::FDA;
  template <class Base> friend class FDABase;

  /// Map atom/residue pair to scalar forces
  std::map<int, std::map<int, Force<real>>> scalar;

  /// Map atom/residue pair to summed vector forces
  std::map<int, std::map<int, Force<Vector>>> summed;

  /// Map atom/residue pair to detailed forces
  std::map<int, std::map<int, DetailedForce>> detailed;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_ */
