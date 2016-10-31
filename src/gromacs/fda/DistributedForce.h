/*
 * DistributedForce.h
 *
 *  Created on: Oct 31, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DISTRIBUTEDFORCE_H_
#define SRC_GROMACS_FDA_DISTRIBUTEDFORCE_H_

#include <vector>
#include "gromacs/math/vectypes.h"

struct DistributedForcePerAtom
{
  /// Atom index
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

struct DistributedForce
{
  /// Indexing table: real atom nr. to index in the pf array; this has length equal to the total nr. of atoms in system
  std::vector<DistributedForcePerAtom> forces;
};

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCE_H_ */
