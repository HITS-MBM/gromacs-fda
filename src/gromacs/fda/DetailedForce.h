/*
 * DetailedForce.h
 *
 *  Created on: Dec 7, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DETAILEDFORCE_H_
#define SRC_GROMACS_FDA_DETAILEDFORCE_H_

#include <array>
#include "PureInteractionType.h"
#include "Vector.h"

namespace fda {

struct DetailedForce
{
  std::array<Vector, static_cast<int>(PureInteractionType::NUMBER)> force;
};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DETAILEDFORCE_H_ */
