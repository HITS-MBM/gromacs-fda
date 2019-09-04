/*
 *  Stress.h
 *
 *  Created on: Jul 27, 2018
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_STRESS_H_
#define SRC_GROMACS_FDA_STRESS_H_

#include <string>
#include <vector>
#include "gromacs/utility/real.h"

namespace fda {

/// Read and write punctual and virial stress files
struct Stress
{
	typedef std::vector<real> StressType;
	typedef std::vector<StressType> StressFrameArrayType;

    Stress(std::string const& filename);

    template <class Comparer>
    bool equal(Stress const& other, Comparer const& comparer) const
    {
    	StressFrameArrayType stress_array1 = this->get_stress();
    	StressFrameArrayType stress_array2 = other.get_stress();

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

    /// Read stress from file
    StressFrameArrayType get_stress() const;

    /// Write stress to file
    void write(std::string const& out_filename, bool out_binary = false) const;

    /// Returns true if the format is binary
    bool get_is_binary() const { return is_binary; }

private:

    std::string filename;

    bool is_binary;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_STRESS_H_ */
