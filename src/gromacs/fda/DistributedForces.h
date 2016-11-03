/*
 * DistributedForces.h
 *
 *  Created on: Oct 31, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_
#define SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_

#include <vector>
#include "ResultType.h"
#include "gromacs/math/vectypes.h"

namespace fda {

struct DistributedForcesSingle
{
  /// Atom or residue index, respectively
  std::vector<int> atom_id;

  /// Scalar distributed forces
  std::vector<real> scalar_force;

  /// Vector distributed forces
  std::vector<rvec> vector_force;

  /// Detailed vector distributed forces
  std::vector<rvec> vector_force_coulomb;
  std::vector<rvec> vector_force_lj;
  std::vector<rvec> vector_force_nb14;
  std::vector<rvec> vector_force_bonds;
  std::vector<rvec> vector_force_angles;
  std::vector<rvec> vector_force_dihedrals;
  std::vector<rvec> vector_force_polar;
};

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
struct DistributedForces
{
  /// Constructor
  DistributedForces(ResultType const& resultType);

  /// Indexing table: real atom nr. to index in the pf array; this has length equal to the total nr. of atoms in system
  std::vector<DistributedForcesSingle> forces;
};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_ */
