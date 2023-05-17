/*
 * GraphTest.cpp
 *
 *  Created on: Oct 17, 2022
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <iostream>
#include <string>
#include <vector>
#include <gtest/gtest.h>
#include "gromacs/gmxana/fda/Graph.h"
#include "testutils/cmdlinetest.h"

using namespace fda_analysis;

//! Test fixture for FDA
class GraphTest : public gmx::test::CommandLineTestBase
{};

TEST_F(GraphTest, empty)
{
    Graph graph;
}

TEST_F(GraphTest, residue)
{
    std::vector<double> forceMatrix{1.0, 1.0, 1.0, 1.0};
    rvec coord_traj[3] = {{1.0, 1.1, 1.2}, {2.0, 2.1, 2.2}, {3.0, 3.1, 3.2}};
    int index[2] = {0, 1};
    int isize = 2;

    Graph graph(forceMatrix, coord_traj, index, isize);

    std::cout << graph << std::endl;
}
