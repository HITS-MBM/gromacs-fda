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

void PunctualStress::write(std::string const& out_filename, bool out_binary) const
{
    auto&& stress_all_frames = get_stress();
    if (out_binary) {
        std::ofstream os(out_filename, std::ifstream::binary);
        if (!os) gmx_fatal(FARGS, "Error opening file.");

        char b = 'b';
        os.write(&b, 1);

        uint nb_atoms = stress_all_frames[0].size();
        os.write(reinterpret_cast<char*>(&nb_atoms), sizeof(uint));
        for (auto&& stress : stress_all_frames) {
            os.write(reinterpret_cast<char*>(&stress[0]), nb_atoms * sizeof(real));
        }
    } else {
        std::ofstream os(out_filename);
        if (!os) gmx_fatal(FARGS, "Error opening file.");

    	os << ResultType::PUNCTUAL_STRESS << std::endl;

        for (auto&& stress : stress_all_frames) {
            if (stress.size() > 0) os << stress[0];
            for (uint i = 1; i < stress.size(); ++i) {
                os << " " << stress[i];
            }
            os << std::endl;
        }
    }
}

PunctualStress::PunctualStressFrameArrayType PunctualStress::get_stress() const
{
    PunctualStressFrameArrayType stress_all_frames;
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
        PunctualStressType stress(syslen);
        for (;;) {
            is.read(reinterpret_cast<char*>(&stress[0]), syslen * sizeof(real));
            stress_all_frames.push_back(stress);
            if (is.tellg() == length) break;
        }
    } else {
        std::ifstream is(filename);
        if (!is) gmx_fatal(FARGS, "Error opening file.");

        std::string token;
        is >> token;
        if (token != "punctual_stress") gmx_fatal(FARGS, "Wrong file type in PunctualStress::get_stress");
        PunctualStressType stress;

        std::string line;
        if (getline(is, line))
        {
            std::stringstream ss(line);
            real value;
            while(ss >> value) {
               stress.push_back(value);
            }
            stress_all_frames.push_back(stress);
        }

        while (getline(is, line))
        {
            std::stringstream ss(line);
            real value;
            while(ss >> value) {
               stress.push_back(value);
            }
            stress_all_frames.push_back(stress);
        }
    }
    return stress_all_frames;
}

} // namespace fda
