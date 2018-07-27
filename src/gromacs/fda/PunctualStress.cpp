/*
 *  PunctualStress.h
 *
 *  Created on: Jul 27, 2018
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "gromacs/utility/fatalerror.h"
#include "PunctualStress.h"
#include "ResultType.h"

namespace fda {

PunctualStress::PunctualStress(std::string const& filename)
 : filename(filename),
   is_binary(false)
{
    std::ifstream file(filename);
    if (!file) gmx_fatal(FARGS, "Error opening file.");
    char first_character;
    file.read(&first_character, 1);
    if (first_character == 'b') is_binary = true;
}

void PunctualStress::write(std::string const& filename, bool binary = false) const
{
	auto&& stress = get_stress();
    if (binary == false) {
    } else {
    }
}

PunctualStress::PunctualStressArray PunctualStress::get_stress() const
{
	PunctualStressArray stress;
    if (this->is_binary) {
    	std::ifstream is(filename, std::ifstream::binary);
		if (!is) gmx_fatal(FARGS, "Error opening file.");

		// get length of file:
		is.seekg (0, is.end);
		int length = is.tellg();
		is.seekg (0, is.beg);

		char first_character;
		is.read(&first_character, 1);
		if (first_character != 'b') gmx_fatal(FARGS, "Wrong file type in PunctualStress::write");

    	uint syslen;
		is.read(reinterpret_cast<char*>(&syslen), sizeof(uint));
		stress.resize(syslen);
		for (;;) {
			is.read(reinterpret_cast<char*>(&i), sizeof(uint));
			is.read(reinterpret_cast<char*>(&nb_interactions_of_i), sizeof(uint));
			for (int m = 0; m != syslen; ++m, ++n) {
				is.read(reinterpret_cast<char*>(&stress), sizeof(real));
				pairwise_forces.push_back(PairwiseForce<ForceType>(i, j, ForceType(force, type)));
			}
		}
    } else {

    }
    return stress;
}

} // namespace fda
