/*
 * PDBTest.cpp
 *
 *  Created on: Mar 18, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "gromacs/gmxana/fda/PDB.h"
#include "gromacs/utility/smalloc.h"
#include "testutils/integrationtests.h"

using namespace fda_analysis;

//! Test fixture for FDA
class PDBTest : public ::gmx::test::IntegrationTestFixture
{};

TEST_F(PDBTest, Basic)
{
	std::vector<int> atomGroup{2, 12, 22};
    PDB pdb(std::string(fileManager_.getInputDataDirectory()) + "/data/glycine_trimer/glycine_trimer.pdb", atomGroup);

    rvec *newCoords;
    snew(newCoords, 30);
    for (int i = 0; i != 30; ++i)
        for (int j = 0; j != 3; ++j)
            newCoords[i][j] = (i*3 + j) * 0.1;

    pdb.updateCoordinates(newCoords);

    for (int i = 0; i != 3; ++i)
        for (int j = 0; j != 3; ++j)
            EXPECT_DOUBLE_EQ(pdb.getCoordinates()[i][j], atomGroup[i]*3 + j);
}
