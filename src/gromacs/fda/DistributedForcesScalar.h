/*
 * DistributedForcesScalar.h
 *
 *  Created on: Nov 4, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DISTRIBUTEDFORCESSCALAR_H_
#define SRC_GROMACS_FDA_DISTRIBUTEDFORCESSCALAR_H_

#include <vector>
#include "gromacs/math/vectypes.h"

namespace fda {

class DistributedForcesScalar
{
public:

  /// Constructor
  DistributedForcesScalar();

  /// Divide all scalar forces by the divisor
  void scalar_real_divide(real divisor);

private:

  /// Atom or residue index, respectively
  std::vector<int> id;

  /// Distributed forces as scalar
  std::vector<real> scalar_forces;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCESSCALAR_H_ */
