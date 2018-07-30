/*
 *  PunctualStress.h
 *
 *  Created on: Jul 27, 2018
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_PUNCTUALSTRESS_H_
#define SRC_GROMACS_FDA_PUNCTUALSTRESS_H_

#include <string>
#include <vector>
#include "gromacs/utility/real.h"

namespace fda {

/**
 * Read pairwise forces from file into arrays and compare.
 */
struct PunctualStress
{
	typedef std::vector<real> PunctualStressType;
	typedef std::vector<PunctualStressType> PunctualStressFrameArrayType;

    PunctualStress(std::string const& filename);

    template <class Comparer>
    bool equal(PunctualStress const& other, Comparer const& comparer) const
    {
    	PunctualStressFrameArrayType stress_array1 = this->get_stress();
    	PunctualStressFrameArrayType stress_array2 = other.get_stress();

        if (stress_array1.size() != stress_array2.size()) return false;

        for (size_t frame = 0; frame != stress_array1.size(); ++frame) {
            auto&& ps1 = stress_array1[frame];
            auto&& ps2 = stress_array1[frame];
            if (ps1.size() != ps2.size()) return false;

            for (size_t i = 0; i != ps1.size(); ++i) {
                if (!comparer(ps1[i], ps2[i])) return false;
            }
        }
        return true;
    }

    /// Write punctual stress to file
    void write(std::string const& out_filename, bool out_binary = false) const;

    /// Returns true if the format is binary
    bool get_is_binary() const { return is_binary; }

private:

    PunctualStressFrameArrayType get_stress() const;

    std::string filename;

    bool is_binary;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_PUNCTUALSTRESS_H_ */
