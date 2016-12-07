/*
 * PairwiseForces.h
 *
 *  Created on: Dec 7, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_PROGRAMS_MDRUN_FDA_TESTS_PAIRWISEFORCES_H_
#define SRC_PROGRAMS_MDRUN_FDA_TESTS_PAIRWISEFORCES_H_

#include <map>
#include <string>
#include <vector>
#include "gromacs/fda/Force.h"
#include "gromacs/fda/Vector.h"
#include "gromacs/utility/real.h"

namespace fda {

/**
 * Read pairwise forces from file into arrays and compare.
 */
template <typename ForceType>
class PairwiseForces
{
public:

	PairwiseForces(std::string const& filename);

	template <class Comparer>
	bool equal(PairwiseForces const& other, Comparer const& comparer) const
	{
		if (all_pairwise_forces.size() != other.all_pairwise_forces.size()) return false;
		for (size_t frame = 0; frame !=  all_pairwise_forces.size(); ++frame) {
            auto pairwise_force = all_pairwise_forces[frame];
            auto other_pairwise_force = other.all_pairwise_forces[frame];
    		if (pairwise_force.size() != other_pairwise_force.size()) return false;
    		for (size_t i = 0; i != pairwise_force.size(); ++i) {
//                auto atom_i = pairwise_force[i];
//                auto other_atom_i = other_pairwise_force[i];
//        		if (atom_i.first != other_atom_i.first) return false;
//        		if (atom_i.second.size() != other_atom_i.second.size()) return false;
    		}
		}
		return true;
	}

private:

	typedef typename std::map<int, std::map<int, ForceType>> PairwiseForceType;

	std::vector<PairwiseForceType> all_pairwise_forces;
};

} // namespace fda

#endif /* SRC_PROGRAMS_MDRUN_FDA_TESTS_PAIRWISEFORCES_H_ */
