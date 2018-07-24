/*
 * PairwiseForces.h
 *
 *  Created on: Dec 7, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_PAIRWISEFORCES_H_
#define SRC_GROMACS_FDA_PAIRWISEFORCES_H_

#include <string>
#include <vector>
#include "gromacs/fda/Force.h"
#include "gromacs/fda/Vector.h"
#include "gromacs/utility/real.h"

namespace fda {

template <typename ForceType>
struct PairwiseForce
{
    PairwiseForce(int i, int j, ForceType force)
     : i(i), j(j), force(force)
    {}

    bool operator == (PairwiseForce const& other) const {
        return i == other.i and j == other.j and force == other.force;
    }

    bool operator != (PairwiseForce const& other) const {
        return !operator == (other);
    }

    template <class Comparer>
    bool equal(PairwiseForce const& other, Comparer const& comparer) const {
        return i == other.i and j == other.j and force.equal(other.force, comparer);
    }

    int i;
    int j;
    ForceType force;
};

/// Output stream
template <typename ForceType>
std::ostream& operator << (std::ostream& os, std::vector<PairwiseForce<ForceType>> const& pairwise_forces)
{
	for (auto&& e : pairwise_forces) {
		os << e.i << " "
		   << e.j << " "
		   << e.force << std::endl;
	}
    return os;
}

/**
 * Read pairwise forces from file into arrays and compare.
 */
template <typename ForceType>
struct PairwiseForces
{
    PairwiseForces(std::string const& filename)
     : filename(filename)
    {}

    template <class Comparer>
    bool equal(PairwiseForces const& other, Comparer const& comparer) const
    {
        std::vector<std::vector<PairwiseForce<ForceType>>> pfl1 = this->get_all_pairwise_forces();
        std::vector<std::vector<PairwiseForce<ForceType>>> pfl2 = other.get_all_pairwise_forces();

        if (pfl1.size() != pfl2.size()) return false;

        // std::zip would be nice; avoid boost::zip as boost was ejected by gromacs
        for (size_t frame = 0; frame != pfl1.size(); ++frame) {
            auto pf1 = pfl1[frame];
            auto pf2 = pfl2[frame];
            if (pf1.size() != pf2.size()) return false;

            for (size_t i = 0; i != pf1.size(); ++i) {
                auto pairwise_force = pf1[i];
                auto other_pairwise_force = pf2[i];
                if (!pairwise_force.equal(other_pairwise_force, comparer)) {
                    std::cout << "i (actual, expected) = " << pairwise_force.i << " " << other_pairwise_force.i << std::endl;
                    std::cout << "j (actual, expected) = " << pairwise_force.j << " " << other_pairwise_force.j << std::endl;
                    std::cout << "force (actual, expected) = " << pairwise_force.force << " " << other_pairwise_force.force << std::endl;
                    return false;
                }
            }
        }
        return true;
    }

    /// Return the number of frames within a pfr-file
    size_t get_number_of_frames() const;

    /// Return pairwise forces of all frames
    std::vector<std::vector<PairwiseForce<ForceType>>> get_all_pairwise_forces(bool sort = false) const;

    /// Sorting the pairwise forces by i, j, and type
    void sort(std::vector<PairwiseForce<ForceType>>& pairwise_forces) const;

    /// Return the maximum index number of the second column and the first frame.
    size_t get_max_index_second_column_first_frame() const;

    /// Parse a file in the scalar format which contains a given number of
    /// atom/residues. Can read a single frame given by the argument frame.
    std::vector<double> get_forcematrix_of_frame(int nbParticles, int frame) const;

    /// Parse a file in the scalar format which contains a given number of
    /// atom/residues. Returns average over all frames.
    std::vector<double> get_averaged_forcematrix(int nbParticles) const;

    /// Write all pairwise forces to file
    void write(std::string const& filename, bool binary = false) const;

    /// Returns true if the format is binary
    bool is_binary() const;

private:

    std::vector<PairwiseForce<ForceType>> get_pairwise_forces(std::ifstream& is) const;
    std::vector<PairwiseForce<ForceType>> get_pairwise_forces_binary(std::ifstream& is) const;

    void write_pairwise_forces(std::ofstream& os, std::vector<PairwiseForce<ForceType>> const& pairwise_forces, int frame) const;
    void write_pairwise_forces_binary(std::ofstream& os, std::vector<PairwiseForce<ForceType>> const& pairwise_forces) const;

    /// Output stream
    friend std::ostream& operator << (std::ostream& os, PairwiseForces const& pf)
    {
        std::vector<std::vector<PairwiseForce<ForceType>>> all_pairwise_forces = pf.get_all_pairwise_forces();
        for (size_t i = 0; i != all_pairwise_forces.size(); ++i) {
            os << "frame " << i << std::endl;
            for (auto const& e : all_pairwise_forces[i]) {
                os << e.i << " "
                   << e.j << " "
                   << e.force << std::endl;
            }
        }
        return os;
    }

    std::string filename;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_PAIRWISEFORCES_H_ */
