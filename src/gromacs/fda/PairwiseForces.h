/*
 * PairwiseForces.h
 *
 *  Created on: Dec 7, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_PAIRWISEFORCES_H_
#define SRC_GROMACS_FDA_PAIRWISEFORCES_H_

#include <map>
#include <string>
#include <vector>
#include <gtest/gtest_prod.h>
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

	/// Default constructor
	PairwiseForces() {}

	/// Constructor reading files
	PairwiseForces(std::string const& filename);

	template <class Comparer>
	bool equal(PairwiseForces const& other, Comparer const& comparer) const
	{
		if (all_pairwise_forces.size() != other.all_pairwise_forces.size()) {
			std::cout << "number of frames       = " << all_pairwise_forces.size() << std::endl;
			std::cout << "number of frames (ref) = " << other.all_pairwise_forces.size() << std::endl;
			return false;
		}
		for (size_t frame = 0; frame !=  all_pairwise_forces.size(); ++frame) {
            auto pairwise_force = all_pairwise_forces[frame];
            auto other_pairwise_force = other.all_pairwise_forces[frame];
    		if (pairwise_force.size() != other_pairwise_force.size()) {
    			std::cout << "frame number = " << frame << std::endl;
    			std::cout << "number of first atoms       = " << pairwise_force.size() << std::endl;
    			std::cout << "number of first atoms (ref) = " << other_pairwise_force.size() << std::endl;
    			return false;
    		}

    		// std::zip would be nice; avoid boost::zip as boost was ejected by gromacs
            for (auto const& atom_i : pairwise_force) {
            	auto other_atom_i = other_pairwise_force.find(atom_i.first);
                if (other_atom_i == other_pairwise_force.end()) return false;
                if (atom_i.second.size() != other_atom_i->second.size()) {
        			std::cout << "frame number = " << frame << std::endl;
        			std::cout << "atom i = " << atom_i.first << std::endl;
        			std::cout << "number of second atoms       = " << atom_i.second.size() << std::endl;
        			std::cout << "number of second atoms (ref) = " << other_atom_i->second.size() << std::endl;
        			return false;
        		}
    		}
		}
		return true;
	}

private:

	FRIEND_TEST(PairwiseForcesTest, DefaultConstructor);
	FRIEND_TEST(PairwiseForcesTest, ReadFile);

	typedef typename std::map<int, std::map<int, ForceType>> PairwiseForceType;

	std::vector<PairwiseForceType> all_pairwise_forces;
};

} // namespace fda

#endif /* SRC_GROMACS_FDA_PAIRWISEFORCES_H_ */
