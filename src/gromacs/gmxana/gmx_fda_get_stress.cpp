/*
 * gmx_fda_get_stress.cpp
 *
 *  Created on: Feb 13, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstddef>
#include <vector>
#include "fda/Graph.h"
#include "fda/Helpers.h"
#include "gmx_ana.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fda/PairwiseForces.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

using namespace fda_analysis;

#define PRINT_DEBUG

int gmx_fda_get_stress(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] calculate the punctual stress by pairwise forces. "
        "If the optional file [TT]-diff[tt] is used "
        "the differences of the pairwise forces will be taken."
    };

    gmx_output_env_t *oenv;

    t_filenm fnm[] = {
        { efPFX, "-i", nullptr, ffREAD },
        { efPFX, "-diff", nullptr, ffOPTRD },
        { efPSX, "-o", nullptr, ffWRITE }
    };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME,
        NFILE, fnm, 0, nullptr, asize(desc), desc, 0, nullptr, &oenv)) return 0;

    // Open pairwise forces file
    fda::PairwiseForces<fda::Force<real>> pairwise_forces(opt2fn("-i", NFILE, fnm));
    fda::PairwiseForces<fda::Force<real>> pairwise_forces_diff(opt2fn("-diff", NFILE, fnm));

    size_t nbFrames = pairwise_forces.get_number_of_frames();
    size_t nbParticles = pairwise_forces.get_max_index_second_column_first_frame() + 1;
    size_t nbParticles2 = nbParticles * nbParticles;

    if (opt2bSet("-diff", NFILE, fnm) and (fn2ftp(opt2fn("-diff", NFILE, fnm)) != fn2ftp(opt2fn("-i", NFILE, fnm))))
        gmx_fatal(FARGS, "Type of the file (-diff) does not match the type of the file (-i).");

    if (opt2bSet("-diff", NFILE, fnm) and pairwise_forces_diff.get_number_of_frames() != nbFrames)
        gmx_fatal(FARGS, "Number of frames is not identical between the two pairwise force files.");

#ifdef PRINT_DEBUG
        std::cout << "pfx filename = " << opt2fn("-i", NFILE, fnm) << std::endl;
        if (opt2bSet("-diff", NFILE, fnm)) std::cout << "pfx-diff filename = " << opt2fn("-diff", NFILE, fnm) << std::endl;
        std::cout << "result filename = " << opt2fn("-o", NFILE, fnm) << std::endl;
        std::cout << "nbFrames = " << nbFrames << std::endl;
        std::cout << "nbParticles = " << nbParticles << std::endl;
#endif

    std::ofstream opsFile(opt2fn("-o", NFILE, fnm));
    if (!opsFile) gmx_fatal(FARGS, "Error opening file %s", opt2fn("-o", NFILE, fnm));
    opsFile << std::scientific << std::setprecision(6)
            << "punctual_stress\n";

    for (size_t frame = 0; frame != nbFrames; ++frame)
    {
        std::vector<double> forceMatrix = pairwise_forces.get_forcematrix_of_frame(nbParticles, frame);
        if (opt2bSet("-diff", NFILE, fnm)) {
            std::vector<double> forceMatrix2 = pairwise_forces_diff.get_forcematrix_of_frame(nbParticles, frame);
            for (size_t i = 0; i < nbParticles2; ++i) forceMatrix[i] -= forceMatrix2[i];
        }

        for (size_t i = 0; i < nbParticles; ++i) {
        real value = 0.0;
            for (size_t j = 0; j < nbParticles; ++j) value += std::abs(forceMatrix[i*nbParticles + j]);
            opsFile << value << " ";
        }
        opsFile << std::endl;
    }

    std::cout << "All done." << std::endl;
    return 0;

}
