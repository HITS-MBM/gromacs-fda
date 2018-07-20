/*
 * PairwiseForces.cpp
 *
 *  Created on: Dec 7, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "gromacs/utility/fatalerror.h"
#include "PairwiseForces.h"

namespace fda {

template <typename ForceType>
size_t PairwiseForces<ForceType>::get_number_of_frames() const
{
    std::ifstream file(filename);
    if (!file) gmx_fatal(FARGS, "Error opening file.");

    std::string line;
    size_t number_of_frames = 0;

    getline(file, line);
    if (line != "pairwise_forces_scalar") gmx_fatal(FARGS, "Wrong file type.");

    while (getline(file, line))
    {
        if (line.find("frame") != std::string::npos) {
            ++number_of_frames;
        }
    }

    return number_of_frames;
}

template <typename ForceType>
std::vector<std::vector<PairwiseForce<ForceType>>> PairwiseForces<ForceType>::get_all_pairwise_forces(bool sort) const
{
    std::vector<std::vector<PairwiseForce<ForceType>>> all_pairwise_forces;
    std::ifstream is(filename);
    std::string token;
    is >> token; // skip file type
    is >> token >> token; // skip first frame
    for (;;)
    {
    	auto&& pairwise_forces = get_pairwise_forces(is);
    	if (pairwise_forces.empty()) break;
    	if (sort) this->sort(pairwise_forces);
    	all_pairwise_forces.push_back(pairwise_forces);
    }
    return all_pairwise_forces;
}

template <typename ForceType>
void PairwiseForces<ForceType>::sort(std::vector<PairwiseForce<ForceType>>& pairwise_forces) const
{
    std::sort(pairwise_forces.begin(), pairwise_forces.end(),
        [](PairwiseForce<ForceType> const& pf1, PairwiseForce<ForceType> const& pf2) {
            return pf1.i < pf2.i or
                  (pf1.i == pf2.i and pf1.j < pf2.j) or
                  (pf1.i == pf2.i and pf1.j == pf2.j and pf1.force.type < pf2.force.type);
        });
}

template <typename ForceType>
size_t PairwiseForces<ForceType>::get_max_index_second_column_first_frame() const
{
    std::ifstream file(filename);
    if (!file) gmx_fatal(FARGS, "Error opening file.");

    std::string line;
    size_t i, j, maxIndex = 0;

    getline(file, line);
    if (line != "pairwise_forces_scalar") gmx_fatal(FARGS, "Wrong file type.");

    getline(file, line); // skip frame line

    while (getline(file, line))
    {
        if (line.find("frame") != std::string::npos) break;

        std::stringstream ss(line);
        ss >> i >> j;

        if (j > maxIndex) maxIndex = j;
    }

    return maxIndex;
}

template <typename ForceType>
std::vector<double> PairwiseForces<ForceType>::get_forcematrix_of_frame(int nbParticles, int frame) const
{
    int nbParticles2 = nbParticles * nbParticles;

    std::vector<double> array(nbParticles2, 0.0);
    std::ifstream file(filename);
    if (!file) gmx_fatal(FARGS, "Error opening file.");

    std::string line;
    bool foundFrame = false;
    int frameNumber;
    std::string tmp;
    int i, j;
    double value;

    getline(file, line);
    if (line != "pairwise_forces_scalar") gmx_fatal(FARGS, "Wrong file type.");

    while (getline(file, line))
    {
        if (line.find("frame") != std::string::npos)
        {
            if (foundFrame) break;
            std::istringstream iss(line);
            iss >> tmp >> frameNumber;
            if (frameNumber != frame) continue;
            foundFrame = true;
            continue;
        }
        if (foundFrame) {
            std::istringstream iss(line);
            iss >> i >> j >> value;
            if (i >= nbParticles or j >= nbParticles)
                gmx_fatal(FARGS, "Index is larger than dimension.");
            array[i*nbParticles+j] = value;
            array[j*nbParticles+i] = value;
        }
    }

    if (!foundFrame) gmx_fatal(FARGS, "Frame not found.");
    return array;
}

template <typename ForceType>
std::vector<double> PairwiseForces<ForceType>::get_averaged_forcematrix(int nbParticles) const
{
    std::ifstream file(filename);
    if (!file) gmx_fatal(FARGS, "Error opening file.");

    int nbParticles2 = nbParticles * nbParticles;
    std::vector<double> forcematrix(nbParticles2, 0.0);

    std::string line;
    int numberOfFrames = 0;
    std::string tmp;
    int i, j;
    double value;

    getline(file, line);
    if (line != "pairwise_forces_scalar") gmx_fatal(FARGS, "Wrong file type.");

    while (getline(file, line))
    {
        if (line.find("frame") != std::string::npos) {
            ++numberOfFrames;
            continue;
        }
        std::istringstream iss(line);
        iss >> i >> j >> value;
        if (i < 0 or i >= nbParticles or j < 0 or j >= nbParticles)
            gmx_fatal(FARGS, "Index error in getAveragedForcematrix.");
        forcematrix[i*nbParticles+j] = value;
    }

    if (!numberOfFrames) gmx_fatal(FARGS, "No frame found.");
    for (auto & elem : forcematrix) elem /= numberOfFrames;

    return forcematrix;
}

template <typename ForceType>
void PairwiseForces<ForceType>::write(std::string const& out_filename, bool out_binary) const
{
    if (this->is_binary() == true and out_binary == false) {
        std::ifstream file(filename, std::ifstream::binary);
        if (!file) gmx_fatal(FARGS, "Error opening file.");
        std::ofstream ofile(out_filename);
        if (!ofile) gmx_fatal(FARGS, "Error opening file.");

    } else if (this->is_binary() == false and out_binary == true) {
        std::ifstream is(filename, std::ifstream::binary);
        if (!is) gmx_fatal(FARGS, "Error opening file.");
        std::ofstream os(out_filename);
        if (!os) gmx_fatal(FARGS, "Error opening file.");

        std::string token;
        for (int i = 0; i != 3; ++i) is >> token;
        for (;;)
        {
        	auto&& pairwise_forces = get_pairwise_forces(is);
        	if (pairwise_forces.empty()) break;
            write_pairwise_forces_binary(os, pairwise_forces);
        }
    } else {
        gmx_fatal(FARGS, "Wrong binary mode in PairwiseForces<ForceType>::write");
    }
}

template <typename ForceType>
bool PairwiseForces<ForceType>::is_binary() const
{
    std::ifstream file(filename);
    if (!file) gmx_fatal(FARGS, "Error opening file.");
    char first_character;
    file.read(&first_character, 1);
    return first_character == 'b';
}

template <typename ForceType>
std::vector<PairwiseForce<ForceType>> PairwiseForces<ForceType>::get_pairwise_forces(std::ifstream& is) const
{
    int i, j;
    ForceType force;
    std::vector<PairwiseForce<ForceType>> pairwise_forces;
    std::string token;
    while (is >> token)
    {
        if (token == "frame") {
            is >> token;
            return pairwise_forces;
        }

        try {
           i = std::stoi(token);
        } catch ( ... ) {
            std::cerr << token << std::endl;
            throw;
        }

        is >> j >> force;
        pairwise_forces.push_back(PairwiseForce<ForceType>(i,j,force));
    }
    return pairwise_forces;
}

template <typename ForceType>
std::vector<PairwiseForce<ForceType>> PairwiseForces<ForceType>::get_pairwise_forces_binary(std::ifstream& is) const
{}

template <typename ForceType>
void PairwiseForces<ForceType>::write_pairwise_forces(std::ofstream& os, std::vector<PairwiseForce<ForceType>> const& pairwise_forces) const
{}

template <typename ForceType>
void PairwiseForces<ForceType>::write_pairwise_forces_binary(std::ofstream& os, std::vector<PairwiseForce<ForceType>> const& pairwise_forces) const
{
	uint nb_interaction = pairwise_forces.size();
    os.write(reinterpret_cast<char*>(&nb_interaction), sizeof(uint));
    auto&& iter_end = pairwise_forces.end();
    for (auto&& iter = pairwise_forces.begin(); iter != iter_end; ++iter) {
        os.write(reinterpret_cast<const char*>(&(iter->i)), sizeof(uint));
        uint nb_interactions_of_i = 0;
        auto&& iter2 = iter + 1;
        for (; iter2 != iter_end; ++iter2) {
        	++nb_interactions_of_i;
        	if (iter2->i != iter->i) break;
        }
        os.write(reinterpret_cast<char*>(&nb_interactions_of_i), sizeof(uint));
        for (; iter != iter2; ++iter) {
            os.write(reinterpret_cast<const char*>(&iter->j), sizeof(uint));
            os.write(reinterpret_cast<const char*>(&iter->force.force), sizeof(real));
            os.write(reinterpret_cast<const char*>(&iter->force.type), sizeof(uint));
        }
    }
}

/// template instantiation
template class PairwiseForces<Force<real>>;
template class PairwiseForces<Force<Vector>>;

} // namespace fda
