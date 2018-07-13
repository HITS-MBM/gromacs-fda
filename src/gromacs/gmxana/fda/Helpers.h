/*
 * Helpers.h
 *
 *  Created on: Feb 12, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef HELPERS_H_
#define HELPERS_H_

#include <iterator>
#include <string>
#include <vector>
#include "FrameType.h"
#include "gromacs/utility/real.h"

namespace fda_analysis {

/**
 * Return stress matrix as vector.
 */
std::vector<real> readStress(std::string const& filename, int& nbFrames,
    int& nbParticles);

/**
 * Return true if string is an integer.
 */
bool isInteger(std::string const& str);

/**
 * Return true if string has extension.
 */
bool hasExtension(std::string const& str, std::string const& extension);

/**
 * Return frameType and value for SINGLE, AVERAGE, and SKIP as argument.
 * For frameType ALL the value is set to 1.
 */
FrameType getFrameTypeAndSkipValue(std::string const& str, int& value);

} // namespace fda_analysis

#endif /* HELPERS_H_ */
