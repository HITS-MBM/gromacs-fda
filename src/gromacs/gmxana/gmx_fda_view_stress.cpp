#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include "fda/FrameType.h"
#include "fda/Graph.h"
#include "fda/Helpers.h"
#include "fda/ParticleType.h"
#include "fda/StressType.h"
#include "gromacs/commandline/pargs.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "gromacs/utility/smalloc.h"
#include "macros.h"
#include "vec.h"
#include "gromacs/fileio/futil.h"
#include "index.h"
#include "xvgr.h"
#include "rmpbc.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "physics.h"
#include "gmx_ana.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/fileio/matio.h"

using namespace fda_analysis;

#define PRINT_DEBUG

int gmx_fda_view_stress(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] illustrates the punctual or von Mises virial stress of FDA "
        "as xpm or pdb-file. "
        "The x-axis of the xpm file represent the particle number and the y-axis the "
        "frame number. "
        "For the pdb file the Bfactor column will be used for the value of "
        "the stress and helps the coloring as a function of the stress magnitude. "
        "The pdb-file can be visualized with an program of your choice. "
        "For both file types, the number of different colors can be set with the "
        "[TT]-nbColors[tt] option. "
    };

    output_env_t oenv;
    static const char* frameString = "average 1";
    static bool convert = false;
    static int nbColors = 10;

    t_pargs pa[] = {
        { "-frame", FALSE, etSTR, {&frameString}, "Specify a single frame number or \"average n\" to take the mean over every n-th frame"
              " or \"skip n\" to take every n-th frame or \"all\" to take all frames (e.g. for movies)" },
        { "-convert", FALSE, etBOOL, {&convert}, "Convert force unit from kJ/mol/nm into pN (only for punctual stress)" },
        { "-nbColors", FALSE, etINT, {&nbColors}, "Number of colors for xpm and pdb files" }
    };

    t_filenm fnm[] = {
		{ efSTR, NULL, NULL, ffREAD },
        { efTPS, NULL, NULL, ffOPTRD },
        { efTRX, "-traj", NULL, ffOPTRD },
        { efNDX, NULL, NULL, ffOPTRD },
        { efVST, "-o", "result", ffWRITE }
    };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_BE_NICE,
        NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

	StressType stressType = PUNCTUAL;
	ParticleType particleType = ATOM;
    if (fn2ftp(opt2fn("-f", NFILE, fnm)) == efPSA) {
        stressType = PUNCTUAL;
        particleType = ATOM;
    } else if (fn2ftp(opt2fn("-f", NFILE, fnm)) == efPSR) {
        stressType = PUNCTUAL;
        particleType = RESIDUE;
    } else if (fn2ftp(opt2fn("-f", NFILE, fnm)) == efVMA) {
        stressType = VIRIAL_VON_MISES;
        particleType = ATOM;
    } else {
    	gmx_fatal(FARGS, "Unkown stress or particle type of stress file.");
    }

    int frameValue = 0;
    FrameType frameType = getFrameTypeAndSkipValue(frameString, frameValue);

    int nbFrames, nbParticles;
    std::vector<real> stressMatrix = readStress(opt2fn("-f", NFILE, fnm), nbFrames, nbParticles);

    // Interactive input of group name for residue model points
    int isize = 0;
    atom_id *index;
    char *grpname;
    if (ftp2fn_null(efNDX, NFILE, fnm)) {
		fprintf(stderr, "\nSelect group for residue model points:\n");
		rd_index(ftp2fn(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);
        if (isize != nbParticles) gmx_fatal(FARGS, "Number of atoms in group does not match number of FDA points.");
    }

	#ifdef PRINT_DEBUG
	    std::cerr << "stress filename = " << opt2fn("-f", NFILE, fnm) << std::endl;
        std::cerr << "result filename = " << opt2fn("-o", NFILE, fnm) << std::endl;
	    std::cerr << "frameType = " << EnumParser<FrameType>()(frameType) << std::endl;
		std::cerr << "frameValue = " << frameValue << std::endl;
	    std::cerr << "stressType = " << EnumParser<StressType>()(stressType) << std::endl;
	    std::cerr << "particleType = " << EnumParser<ParticleType>()(particleType) << std::endl;
	    std::cerr << "convert = " << convert << std::endl;
    #endif

    if (fn2ftp(opt2fn("-o", NFILE, fnm)) == efPDB and !ftp2fn_null(efTPS, NFILE, fnm)) gmx_fatal(FARGS, "Input structure is missing.");
    if (stressType != PUNCTUAL and convert) gmx_fatal(FARGS, "Option -convert makes only sense for punctual stress.");

	for (auto & elem : stressMatrix) elem = std::abs(elem);

	// Convert from kJ/mol/nm into pN
	if (convert) for (auto & elem : stressMatrix) elem *= 1.66;

	// stressMatrix2 will not used for frameType == SINGLE
	real **stressMatrix2 = NULL;
    int nbFramesAfterSkipping = 1;
	if (frameType != SINGLE) {
        snew(stressMatrix2, nbParticles);
        nbFramesAfterSkipping = ceil(static_cast<real>(nbFrames) / frameValue);
        for (int i = 0; i < nbParticles; ++i) {
            snew(stressMatrix2[i], nbFramesAfterSkipping);
            if (frameType == AVERAGE) {
                for (int j = 0, js = 0; j < nbFramesAfterSkipping; ++j) {
                    real value = 0.0;
                    for (int k = 0; k < frameValue and js < nbFrames; ++k, ++js) value += stressMatrix[js*nbParticles + i];
                    value /= frameValue;
                    stressMatrix2[i][j] = value;
                }
            } else {
                for (int j = 0, js = 0; j < nbFramesAfterSkipping; ++j, js += frameValue) {
                    stressMatrix2[i][j] = stressMatrix[js*nbParticles + i];
                }
            }
        }
	}

	std::string title;
	if (stressType == PUNCTUAL) title += "Punctual ";
	else if (stressType == VIRIAL) title += "Virial ";
	else if (stressType == VIRIAL_VON_MISES) title += "Von Mises virial ";
	title += "stress";
	if (particleType == ATOM) title += " over atoms";
	else if (particleType == RESIDUE) title += " over residues";

	if (fn2ftp(opt2fn("-o", NFILE, fnm)) == efPDB) {
		if (frameType == SINGLE) {
			char buf[256];
			rvec *xp;
			t_topology top;
			int ePBC;
			matrix box;

			read_tps_conf(ftp2fn(efTPS, NFILE, fnm), buf, &top, &ePBC, &xp, NULL, box, TRUE);

			FILE *fp = gmx_ffopen(opt2fn("-o", NFILE, fnm), "w");

			for (int i = 0; i < nbParticles; ++i)
			    top.atoms.pdbinfo[i].bfac = stressMatrix[frameValue*nbParticles + i];

			write_pdbfile(fp, title.c_str(), &top.atoms, xp, ePBC, box, ' ', 0, NULL, TRUE);
			gmx_ffclose(fp);
		} else {
			gmx_fatal(FARGS, "Not implemented yet.");
//			if (!read_first_frame(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&fr,TRX_READ_X))
//		    gmx_fatal(FARGS,"Could not read a frame from %s",ftp2fn(efTRX,NFILE,fnm));
		}
	} else if (fn2ftp(opt2fn("-o", NFILE, fnm)) == efXPM) {
		FILE *out = gmx_ffopen(opt2fn("-o", NFILE, fnm), "w");
	    t_rgb rlo = {1, 1, 1}, rhi = {0, 0, 0};
	    write_xpm(out, 0, title.c_str(), "", "Particle", "Frame", nbParticles, nbFramesAfterSkipping,
		    NULL, NULL, stressMatrix2, *std::min_element(stressMatrix.begin(), stressMatrix.end()),
			*std::max_element(stressMatrix.begin(), stressMatrix.end()), rlo, rhi, &nbColors);
        gmx_ffclose(out);
	} else
	    gmx_fatal(FARGS, "Missing output filename -opdb or -oxpm.");

    std::cout << "All done" << std::endl;
    return 0;

}
