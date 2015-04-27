#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include "fda/BoostGraph.h"
#include "fda/Helpers.h"
#include "fda/PDB.h"
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
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/utility/cstringutil.h"

using namespace fda_analysis;

#define PRINT_DEBUG

int gmx_fda_shortest_path(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] calculate the k-shortest paths between a source node "
    	"and a destination node of a FDA force network. The graph will "
        "printed in the PDB-format, which allow an easy visualization with "
    	"an program of your choice. "
    	"Each path will be determined and segment names will be assign to each "
        "of them, thus coloring them by segment id will help the analysis "
        "(32 different colors). The Bfactor column will be used for the value of "
        "the force and helps the coloring as a function of the force magnitude. "
        "The CONNECT header will be used to create bonds between nodes. "
    };

    output_env_t oenv;
    static const char* frameString = "average 1";
    static int source = 0;
    static int dest = 0;
    static int numberOfShortestPaths = 1;
    static bool convert = false;

    t_pargs pa[] = {
        { "-frame", FALSE, etSTR, {&frameString}, "Specify a single frame number or \"average n\" to take the mean over every n-th frame"
                " or \"skip n\" to take every n-th frame or \"all\" to take all frames (e.g. for movies)" },
        { "-source", FALSE, etINT, {&source}, "Source node of the path" },
        { "-dest", FALSE, etINT, {&dest}, "Destination point of the path" },
        { "-nk", FALSE, etINT, {&numberOfShortestPaths}, "Number of shortest paths" },
        { "-convert", FALSE, etBOOL, {&convert}, "Convert force unit from kJ/mol/nm into pN" }
    };

    t_filenm fnm[] = {
        { efPFX, "-ipf", NULL, ffREAD },
        { efPFX, "-ipf-diff", NULL, ffOPTRD },
        { efTPS, NULL, NULL, ffREAD },
        { efTRX, "-traj", NULL, ffOPTRD },
        { efNDX, NULL, NULL, ffOPTRD },
        { efPDB, "-o", "result", ffWRITE }
    };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_BE_NICE,
        NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv)) return 0;

    if (opt2bSet("-ipf-diff", NFILE, fnm) and (fn2ftp(opt2fn("-ipf-diff", NFILE, fnm)) != fn2ftp(opt2fn("-ipf", NFILE, fnm))))
        gmx_fatal(FARGS, "Type of the file (-ipf-diff) does not match the type of the file (-ipf).");

    // Get number of particles
    int nbParticles = getMaxIndexSecondColumnFirstFrame(opt2fn("-ipf", NFILE, fnm)) + 1;
    int nbParticles2 = nbParticles * nbParticles;

    // Interactive input of group name for residue model points
    int isize = 0;
    atom_id *index;
    char *grpname;
    fprintf(stderr, "\nSelect group for residue model points:\n");
    rd_index(ftp2fn(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);
    if (isize != nbParticles) gmx_fatal(FARGS, "Number of atoms in group does not match number of FDA points.");

    int frameValue;
    FrameType frameType = getFrameTypeAndSkipValue(frameString, frameValue);

	#ifdef PRINT_DEBUG
		std::cerr << "frameType = " << EnumParser<FrameType>()(frameType) << std::endl;
		std::cerr << "frameValue = " << frameValue << std::endl;
		std::cerr << "Number of particles (np) = " << nbParticles << std::endl;
		std::cerr << "source = " << source << std::endl;
		std::cerr << "dest = " << dest << std::endl;
		std::cerr << "Number of shortest paths (nk) = " << numberOfShortestPaths << std::endl;
		std::cerr << "convert = " << convert << std::endl;
		std::cerr << "pfx filename = " << opt2fn("-ipf", NFILE, fnm) << std::endl;
		if (opt2bSet("-ipf-diff", NFILE, fnm)) std::cerr << "pfx-diff filename = " << opt2fn("-ipf-diff", NFILE, fnm) << std::endl;
	    std::cerr << "structure filename = " << opt2fn("-s", NFILE, fnm) << std::endl;
		std::cerr << "result filename = " << opt2fn("-o", NFILE, fnm) << std::endl;
		std::cerr << "trajectory filename = " << opt2fn("-traj", NFILE, fnm) << std::endl;
    #endif

    // Read input structure coordinates
    char buf[256];
    rvec *coord;
    t_topology top;
    int ePBC;
    matrix box;
    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), buf, &top, &ePBC, &coord, NULL, box, TRUE);

    PDB pdb(opt2fn("-s", NFILE, fnm), std::vector<int>(index, index + isize));

    if (frameType == ALL)
    {
        // Read trajectory coordinates
        t_trxstatus *status;
        real time;
        rvec *coord_traj;
        matrix box;
        read_first_x(oenv, &status, opt2fn("-traj", NFILE, fnm), &time, &coord_traj, box);

        int nbFrames = getNumberOfFrames(opt2fn("-ipf", NFILE, fnm));
        for (int frame = 0; frame != nbFrames; frame += frameValue)
        {
    		std::vector<double> forceMatrix, forceMatrix2;
			forceMatrix = parseScalarFileFormat(opt2fn("-ipf", NFILE, fnm), nbParticles, frame);
			if (opt2bSet("-ipf-diff", NFILE, fnm)) forceMatrix2 = parseScalarFileFormat(opt2fn("-ipf-diff", NFILE, fnm), nbParticles, frame);

    		if (opt2bSet("-ipf-diff", NFILE, fnm)) for (int i = 0; i < nbParticles2; ++i) forceMatrix[i] -= forceMatrix2[i];
    		for (auto & f : forceMatrix) f = std::abs(f);

    		// Convert from kJ/mol/nm into pN
    		if (convert) for (auto & f : forceMatrix) f *= 1.66;

    		BoostGraph graph(forceMatrix);
    		BoostGraph::PathList shortestPaths = graph.findKShortestPaths(source, dest, numberOfShortestPaths);

            read_next_x(oenv, status, &time, coord_traj, box);
    		pdb.updateCoordinates(coord_traj);
    		pdb.writePaths(opt2fn("-o", NFILE, fnm), shortestPaths, forceMatrix, frame);
        }

        close_trj(status);

    } else {

		std::vector<double> forceMatrix, forceMatrix2;

		if (frameType == SINGLE) {
			int frame = atoi(frameString);
			forceMatrix = parseScalarFileFormat(opt2fn("-ipf", NFILE, fnm), nbParticles, frame);
			if (opt2bSet("-ipf-diff", NFILE, fnm)) forceMatrix2 = parseScalarFileFormat(opt2fn("-ipf-diff", NFILE, fnm), nbParticles, frame);
		} else if (frameType == AVERAGE) {
			forceMatrix = getAveragedForcematrix(opt2fn("-ipf", NFILE, fnm), nbParticles);
			if (opt2bSet("-ipf-diff", NFILE, fnm)) forceMatrix2 = getAveragedForcematrix(opt2fn("-ipf", NFILE, fnm), nbParticles);
		} else {
			gmx_fatal(FARGS, "Unknown frame type: %s", EnumParser<FrameType>()(frameType).c_str());
		}

		if (opt2bSet("-ipf-diff", NFILE, fnm)) for (int i = 0; i < nbParticles2; ++i) forceMatrix[i] -= forceMatrix2[i];
		for (auto & f : forceMatrix) f = std::abs(f);

		// Convert from kJ/mol/nm into pN
		if (convert) for (auto & f : forceMatrix) f *= 1.66;

		BoostGraph graph(forceMatrix);
		BoostGraph::PathList shortestPaths = graph.findKShortestPaths(source, dest, numberOfShortestPaths);

		pdb.writePaths(opt2fn("-o", NFILE, fnm), shortestPaths, forceMatrix, false);
    }

    std::cout << "All done" << std::endl;
    return 0;

}
