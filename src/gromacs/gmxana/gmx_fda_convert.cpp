/*
 * fda_convert.cpp
 *
 *  Created on: Jul 16, 2018
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <fstream>
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

using namespace fda;

int main(int argc, char *argv[])
{
    std::cout << "FDA convert" << std::endl;

    const char *desc[] = {
        "[THISMODULE] converts pairwise forces and punctual stress from text- into binary-format and vice versa."
    };

    gmx_output_env_t *oenv;

    t_pargs pa[] = {};

    t_filenm fnm[] = {
        { efPFX, "-i", NULL, ffREAD },
        { efPFX, "-o", NULL, ffWRITE }
    };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME,
        NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv)) return 0;

    // Open pairwise forces file
    fda::PairwiseForces<fda::Force<real>> pairwise_forces(opt2fn("-i", NFILE, fnm));

#ifdef PRINT_DEBUG
        std::cout << "input file = " << opt2fn("-i", NFILE, fnm) << std::endl;
        std::cout << "output file = " << opt2fn("-o", NFILE, fnm) << std::endl;
#endif

    std::ofstream opsFile(opt2fn("-o", NFILE, fnm));
    if (!opsFile) gmx_fatal(FARGS, "Error opening file", opt2fn("-o", NFILE, fnm));
    opsFile << std::scientific << std::setprecision(6)
            << "punctual_stress\n";

    std::cout << "All done." << std::endl;
    return 0;
}
