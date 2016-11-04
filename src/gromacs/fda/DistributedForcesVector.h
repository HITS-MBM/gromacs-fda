/*
 * DistributedForcesVector.h
 *
 *  Created on: Nov 4, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DISTRIBUTEDFORCESVECTOR_H_
#define SRC_GROMACS_FDA_DISTRIBUTEDFORCESVECTOR_H_

#include <vector>
#include "gromacs/math/vectypes.h"

namespace fda {

class DistributedForcesVector
{
public:

  /// Constructor
  DistributedForcesVector();

private:

  /// Atom or residue index, respectively
  std::vector<int> id;

  /// Distributed forces as vector
  std::vector<rvec> vector_forces;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCESVECTOR_H_ */
