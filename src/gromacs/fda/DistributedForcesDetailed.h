/*
 * DistributedForcesDetailed.h
 *
 *  Created on: Nov 4, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DISTRIBUTEDFORCESDETAILED_H_
#define SRC_GROMACS_FDA_DISTRIBUTEDFORCESDETAILED_H_

#include <array>
#include <map>
#include <vector>
#include "InteractionType.h"
#include "Vector.h"

namespace fda {

class DistributedForcesDetailed
{
public:

  /// Constructor
  DistributedForcesDetailed();

  /**
   * Looks for the interaction in the array corresponding to type; if found,
   * the force is added to the existing value; if not, a new interaction is
   * appended to the array.
   */
  void add_force(int j, int type, rvec force);

private:

  /// Atom or residue index, respectively
  std::vector<int> id;

  /// Distributed forces for each interaction type
  std::array<std::map<int, Vector>, InteractionType::num> interactions;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCESDETAILED_H_ */
