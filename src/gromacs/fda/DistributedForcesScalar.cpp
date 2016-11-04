/*
 * DistributedForcesScalar.cpp
 *
 *  Created on: Nov 4, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "DistributedForcesScalar.h"

using namespace fda;

DistributedForcesScalar::DistributedForcesScalar()
{

}

void DistributedForcesScalar::scalar_real_divide(real divisor)
{
  real inv = 1 / divisor;
  for (auto& f : scalar_forces) f *= inv;
}
