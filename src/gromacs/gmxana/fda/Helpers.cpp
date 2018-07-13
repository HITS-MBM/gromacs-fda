/*
 * Helpers.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "Helpers.h"
#include "gromacs/utility/fatalerror.h"
#include <cctype>
#include <iostream>
#include <fstream>
#include <sstream>

namespace fda_analysis {

std::vector<real> readStress(std::string const& filename, int& nbFrames,
    int& nbParticles)
{
    std::ifstream file(filename);
    if (!file) gmx_fatal(FARGS, "Error opening file %s.", filename.c_str());

    std::string line;
    std::vector<real> stress;
    nbFrames = 0;
    nbParticles = 0;

    getline(file, line);
    if (line != "punctual_stress" and line != "virial_stress_von_mises") gmx_fatal(FARGS, "Wrong file type.");

    if (getline(file, line))
    {
		std::stringstream ss(line);
		real value;
		while(ss >> value) {
		   stress.push_back(value);
			++nbParticles;
		}
		++nbFrames;
    }

    while (getline(file, line))
    {
		std::stringstream ss(line);
		real value;
		while(ss >> value) {
		   stress.push_back(value);
		}
		++nbFrames;
    }

    return stress;
}

bool isInteger(std::string const& str)
{
	for (auto c : str) if (!std::isdigit(c)) return false;
	return true;
}

bool hasExtension(std::string const& str, std::string const& extension)
{
	return str.size() > extension.size() and str.compare(str.size() - extension.size(), extension.size(), extension) == 0;
}

FrameType getFrameTypeAndSkipValue(std::string const& frameString, int& value)
{
	FrameType frameType;
	value = 1;

	std::istringstream iss(frameString);
	std::istream_iterator<std::string> iterCur(iss), iterEnd;
	if (iterCur == iterEnd) gmx_fatal(FARGS, "Error in frame option.");

	std::string firstToken = *iterCur;
    ++iterCur;
	if (iterCur == iterEnd) {
		if (isInteger(firstToken)) {
			frameType = SINGLE;
			value = std::atoi(firstToken.c_str());
		}
		else frameType = EnumParser<FrameType>()(firstToken);
	} else {
	    frameType = EnumParser<FrameType>()(firstToken);
		value = std::atoi(iterCur->c_str());

		// Abort if more than two tokens are found
	    ++iterCur;
		if (iterCur != iterEnd) gmx_fatal(FARGS, "Error in frame option.");
	}

	if (frameType == ALL and value != 1) gmx_fatal(FARGS, "Frame type \"all\" does not expect a number.");

	return frameType;
}

} // namespace fda_analysis
