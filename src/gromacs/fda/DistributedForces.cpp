/*
 * DistributedForces.cpp
 *
 *  Created on: Nov 3, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <DistributedForces.h>

using namespace fda;

DistributedForces::DistributedForces(ResultType const& resultType)
{

}

void DistributedForces::scalar_real_divide(real divisor)
{
  for (auto& f : scalar_forces) {
    f.scalar_real_divide(divisor);
  }
}

void DistributedForces::summed_merge_to_scalar(const rvec *x, int Vector2Scalar)
{
  if (scalar_forces.size() != vector_forces.size())
    gmx_fatal(FARGS, "Mismatch of PF atom id: summed %d vs scalar %d\n", vector_forces.size(), scalar_forces.size());

  for (int i = 0; i != scalar_forces.size(); ++i) {
    scalar_add(scalar_forces, scalar_forces[i].id, scalar_forces[i].type,
      vector2signedscalar(scalar_forces[i].force, x[scalar_forces->nr], x[scalar_forces[i].id], Vector2Scalar));
  }

// Old code:
//  t_pf_interaction_array_summed *ia;
//  t_pf_interaction_summed *i;
//  int j;
//
//  if (src->nr != dst->nr)
//    gmx_fatal(FARGS, "Mismatch of PF atom id: summed %d vs scalar %d\n", src->nr, dst->nr);
//  ia = &(src->interactions);
//  for (j = 0; j < ia->len; j++) {
//    i = &(ia->array[j]);
//    pf_atom_scalar_add(dst, i->jjnr, i->type, pf_vector2signedscalar(i->force, x[src->nr], x[i->jjnr], Vector2Scalar));
//  }
}
