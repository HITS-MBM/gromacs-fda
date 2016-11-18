/*
 * DistributedForcesDetailed.cpp
 *
 *  Created on: Nov 4, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "DistributedForcesDetailed.h"

using namespace fda;

DistributedForcesDetailed::DistributedForcesDetailed()
{

}

void DistributedForcesDetailed::add_force(int j, int type, rvec force)
{
  //interactions[type][j] += force;

  // old:
//  t_pf_interaction_array_detailed *ia;
//  t_pf_interaction_detailed *p;
//
//  ia = pf_interaction_array_by_type(&atom->interactions, type);
//  p = pf_lookup_interaction_detailed(ia, jjnr);
//  if (p == NULL) {
//    pf_interaction_array_detailed_append(ia, jjnr, force);
//  } else {
//    rvec_inc(p->force, force);
//  }
}
