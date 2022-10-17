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

TEST_F(GraphTest, text)
{
    Graph graph;
    //Graph graph(forceMatrix, coord_traj, index, isize);
}
