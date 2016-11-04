/*
 * DistributedForcesDetailed.h
 *
 *  Created on: Nov 4, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DISTRIBUTEDFORCESDETAILED_H_
#define SRC_GROMACS_FDA_DISTRIBUTEDFORCESDETAILED_H_

#include <vector>
#include "gromacs/math/vectypes.h"

namespace fda {

class DistributedForcesDetailed
{
public:

  /// Constructor
  DistributedForcesDetailed();

private:

  /// Atom or residue index, respectively
  std::vector<int> id;

  /// Distributed forces as vector for each interaction
  std::vector<rvec> vector_forces_coulomb;
  std::vector<rvec> vector_forces_lj;
  std::vector<rvec> vector_forces_nb14;
  std::vector<rvec> vector_forces_bonds;
  std::vector<rvec> vector_forces_angles;
  std::vector<rvec> vector_forces_dihedrals;
  std::vector<rvec> vector_forces_polar;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCESDETAILED_H_ */
