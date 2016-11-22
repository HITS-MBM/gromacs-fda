/*
 * DistributedForces.cpp
 *
 *  Created on: Nov 3, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "DistributedForces.h"
#include "gromacs/utility/fatalerror.h"
#include "Utilities.h"

using namespace fda;

DistributedForces::DistributedForces(ResultType result_type, OnePair one_pair, int syslen)
 : result_type(result_type),
   one_pair(one_pair),
   syslen(syslen)
{
  // Allocates and fills with -1 the indexing table real atom number to pf number
  if (result_type != ResultType::NO) sys2pf.resize(syslen, -1);
}

void DistributedForces::scalar_real_divide(real divisor)
{
  real inv = 1 / divisor;
  for (auto& si : scalar)
    for (auto& sj : si.second) sj.second *= inv;
}

void DistributedForces::summed_merge_to_scalar(const rvec *x, Vector2Scalar v2s)
{
  for (auto& si : summed) {
	int i = si.first;
    for (auto& sj : si.second) {
      int j = sj.first;
      scalar[i][j] += vector2signedscalar(sj.second.get_rvec(), x[i], x[j], v2s);
    }
  }
}
