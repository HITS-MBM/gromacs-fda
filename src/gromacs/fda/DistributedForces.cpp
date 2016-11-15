/*
 * DistributedForces.cpp
 *
 *  Created on: Nov 3, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <DistributedForces.h>

using namespace fda;

DistributedForces::DistributedForces(ForceType force_type, ResultType result_type, int syslen)
 : force_type(force_type),
   result_type(result_type),
   syslen(syslen)
{
  std::cout << "Allocating space for pairwise forces of " << force_type << ": " << syslen << std::endl;
  switch(fda_settings.one_pair) {
	case OnePair::DETAILED:
	  snew(atoms->detailed, atoms->len);
	  for (i = 0; i < syslen; i++)
		if (atoms->sys2pf[i] != -1)
		  atoms->detailed[atoms->sys2pf[i]].nr = i;
	  break;
	case OnePair::SUMMED:
	  snew(atoms->summed, atoms->len);
	  for (i = 0; i < syslen; i++)
		if (atoms->sys2pf[i] != -1)
		  atoms->summed[atoms->sys2pf[i]].nr = i;
	  break;
	default:
	   /* this is the first use of OnePair; calling gmx_fatal on further uses will not happen anymore, as the test here is supposed to catch such cases */
	   gmx_fatal(FARGS, "Unknown value for pf OnePair: 0x%x.\n", OnePair);
	   break;
  }

  // Allocates and fills with -1 the indexing table real atom number to pf number
  if (result_type != ResultType::NO) sys2pf.resize(syslen, -1);
}

void DistributedForces::scalar_real_divide(real divisor)
{
  for (auto& f : scalar) f.scalar_real_divide(divisor);
}

void DistributedForces::summed_merge_to_scalar(const rvec *x, int Vector2Scalar)
{
  if (scalar.size() != vector.size())
    gmx_fatal(FARGS, "Mismatch of PF atom id: summed %d vs scalar %d\n", vector.size(), scalar.size());

  for (int i = 0; i != scalar.size(); ++i) {
    scalar_add(scalar, scalar[i].id, scalar[i].type,
      vector2signedscalar(scalar[i].force, x[scalar->nr], x[scalar[i].id], Vector2Scalar));
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
