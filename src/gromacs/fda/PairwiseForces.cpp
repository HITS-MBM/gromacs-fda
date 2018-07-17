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
    int i, j;
    ForceType force;
    std::vector<PairwiseForce<ForceType>> pairwise_forces;
    std::vector<std::vector<PairwiseForce<ForceType>>> all_pairwise_forces;
    std::ifstream is(filename);
    std::string token;
    is >> token;
    while (is >> token)
    {
        if (token == "frame") {
            size_t frame;
            is >> frame;
            if (frame) {
                all_pairwise_forces.push_back(pairwise_forces);
                pairwise_forces.clear();
            }
            if (frame != all_pairwise_forces.size()) {
                std::cout << "frame = " << frame << std::endl;
                std::cout << "all_pairwise_forces.size() = " << all_pairwise_forces.size() << std::endl;
                throw std::runtime_error("frame numbers are not consecutively");
            }
            continue;
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
    all_pairwise_forces.push_back(pairwise_forces);

    if (sort) {
        for (auto & e : all_pairwise_forces) this->sort(e);
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
void PairwiseForces<ForceType>::write(std::string const& ofilename, bool binary = false) const
{
    if (binary) {
    	std::ifstream file(filename, std::ifstream::binary);
    	if (!file) gmx_fatal(FARGS, "Error opening file.");
    	std::ofstream ofile(ofilename);
    	if (!ofile) gmx_fatal(FARGS, "Error opening file.");

    	int i, j;
    	ForceType force;
        while (getline(file, line))
        {
            std::istringstream iss(line);
            iss >> i >> j >> force;


        }
    } else {

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

/// template instantiation
template class PairwiseForces<Force<real>>;
template class PairwiseForces<Force<Vector>>;

} // namespace fda
