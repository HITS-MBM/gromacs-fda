/*
 * PairwiseForces.cpp
 *
 *  Created on: Dec 7, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <fstream>
#include <iostream>
#include <stdexcept>
#include "PairwiseForces.h"

namespace fda {

template <typename ForceType>
PairwiseForces<ForceType>::PairwiseForces(std::string const& filename)
{
	int i, j;
	ForceType force;
	PairwiseForceType pairwise_forces;
	std::ifstream is(filename);
    std::string token;
    while (is >> token)
    {
		if (token == "frame") {
			int frame_number;
			is >> frame_number;
			if (frame_number != all_pairwise_forces.size()) {
				std::cout << "frame_number = " << frame_number << std::endl;
				std::cout << "all_pairwise_forces.size() = " << all_pairwise_forces.size() << std::endl;
				throw std::runtime_error("frame numbers are not consecutively");
			}
			all_pairwise_forces.push_back(pairwise_forces);
			pairwise_forces.clear();
			continue;
		}

		is >> i >> j >> force;
		pairwise_forces[i][j] = force;
    }
}

/// template instantiation
template class PairwiseForces<Force<real>>;
template class PairwiseForces<Force<Vector>>;

} // namespace fda
