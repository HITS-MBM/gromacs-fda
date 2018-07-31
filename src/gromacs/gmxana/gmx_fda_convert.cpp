/*
 * gmx_fda_convert.cpp
 *
 *  Created on: Jul 16, 2018
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <fstream>
#include "gmx_ana.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fda/PairwiseForces.h"
#include <gromacs/fda/Stress.h>
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

#define PRINT_DEBUG 1

int gmx_fda_convert(int argc, char *argv[])
{
    std::cout << "FDA convert" << std::endl;

    const char *desc[] = {
        "[THISMODULE] converts pairwise forces, punctual, and virial stress files"
    	"from text- into binary-format and vice versa."
        "If the input is binary format the output will be text-based and vice versa."
    };

    gmx_output_env_t *oenv;

    t_pargs pa[] = {};

    t_filenm fnm[] = {
        { efFDACON, "-i", NULL, ffREAD },
        { efFDACON, "-o", NULL, ffWRITE }
    };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME,
        NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv)) return 0;

#ifdef PRINT_DEBUG
        std::cout << "input file = " << opt2fn("-i", NFILE, fnm) << std::endl;
        std::cout << "output file = " << opt2fn("-o", NFILE, fnm) << std::endl;
#endif

    if (fn2ftp(opt2fn("-i", NFILE, fnm)) != fn2ftp(opt2fn("-o", NFILE, fnm)))
        gmx_fatal(FARGS, "Input and output file type must be identical.");

    if (fn2ftp(opt2fn("-i", NFILE, fnm)) == efPFA or fn2ftp(opt2fn("-i", NFILE, fnm)) == efPFR) {
        fda::PairwiseForces<fda::Force<real>> pairwise_forces(opt2fn("-i", NFILE, fnm));
        pairwise_forces.write(opt2fn("-o", NFILE, fnm), !pairwise_forces.get_is_binary());
    } else if (fn2ftp(opt2fn("-i", NFILE, fnm)) == efPSA or fn2ftp(opt2fn("-i", NFILE, fnm)) == efPSR
        or fn2ftp(opt2fn("-i", NFILE, fnm)) == efVSA or fn2ftp(opt2fn("-i", NFILE, fnm)) == efVMA) {
        fda::Stress stress(opt2fn("-i", NFILE, fnm));
        stress.write(opt2fn("-o", NFILE, fnm), !stress.get_is_binary());
    }

    std::cout << "All done." << std::endl;
    return 0;
}
