/*
 * FDAGraphTest.cpp
 *
 *  Created on: Feb 4, 2015
 *      Author: Bernd Doser, HITS gGmbH
 */

#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/fileio/path.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "programs/mdrun/mdrun_main.h"
#include "testutils/cmdlinetest.h"
#include "testutils/integrationtests.h"
#include "testutils/TextSplitter.h"
#include "testutils/LogicallyErrorComparer.h"
#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct TestDataStructure
{
	TestDataStructure(
	    std::string const& testDirectory,
	    std::vector<std::string> const& cmdline,
	    std::string const& groupname,
        std::string const& result,
	    std::string const& reference
	)
     : testDirectory(testDirectory),
       cmdline(cmdline),
       groupname(groupname),
       result(result),
	   reference(reference)
	{}

    std::string testDirectory;
    std::vector<std::string> cmdline;
    std::string groupname;
    std::string result;
    std::string reference;
};

} // namespace anonymous

//! Test fixture for FDA
class FDAGraphTest : public ::testing::WithParamInterface<TestDataStructure>,
                     public ::gmx::test::IntegrationTestFixture
{};

//! Test body for FDA
TEST_P(FDAGraphTest, Basic)
{
    std::string cwd = gmx::Path::getWorkingDirectory();
    std::string dataPath = std::string(fileManager_.getInputDataDirectory()) + "/data";
    std::string testPath = fileManager_.getTemporaryFilePath("/" + GetParam().testDirectory);

    std::string cmd = "mkdir -p " + testPath;
    ASSERT_FALSE(system(cmd.c_str()));

    cmd = "cp -r " + dataPath + "/" + GetParam().testDirectory + "/* " + testPath;
    ASSERT_FALSE(system(cmd.c_str()));

    gmx_chdir(testPath.c_str());

    ::gmx::test::CommandLine caller;
    caller.append("gmx_fda fda_graph");
    for (std::vector<std::string>::const_iterator iterCur(GetParam().cmdline.begin()), iterNext(GetParam().cmdline.begin() + 1),
        iterEnd(GetParam().cmdline.end()); iterCur != iterEnd; ++iterCur, ++iterNext)
    {
        if (iterNext == iterEnd or iterNext->substr(0,1) == "-") caller.append(*iterCur);
        else {
            caller.addOption(iterCur->c_str(), iterNext->c_str());
            ++iterCur, ++iterNext;
        }
    }
    caller.addOption("-o", GetParam().result);

    std::cout << caller.toString() << std::endl;

    if (!GetParam().groupname.empty()) redirectStringToStdin((GetParam().groupname + "\n").c_str());

    ASSERT_FALSE(gmx_fda_graph(caller.argc(), caller.argv()));

    const double error_factor = 1.0e4;
    const bool weight_by_magnitude = false;
    const bool ignore_sign = true;

    LogicallyEqualComparer<weight_by_magnitude,ignore_sign> comparer(error_factor);

    // compare atom pairs
    EXPECT_TRUE((compare(TextSplitter(GetParam().reference), TextSplitter(GetParam().result), comparer)));

	gmx_chdir(cwd.c_str());
}

INSTANTIATE_TEST_CASE_P(AllFDAGraphTests, FDAGraphTest, ::testing::Values(
	TestDataStructure(
        "maxime_all_prot",
        {"-ipf", "cap0_all_prot.pfr", "-ipf-diff", "cap1_all_prot.pfr", "-s", "1G6N.pdb", "-n", "index.ndx", "-frame", "0", "-t", "100", "-min", "2", "-convert"},
		"C-alpha",
        "result.pdb",
        "ref.pdb"
    ),
    TestDataStructure(
        "maxime_all_prot",
        {"-ipf", "cap0_all_prot.pfr", "-ipf-diff", "cap1_all_prot.pfr", "-s", "1G6N.pdb", "-n", "index.ndx", "-frame", "0", "-t", "20", "-min", "2", "-convert"},
        "C-alpha",
        "result.dmc",
        "ref.dmc"
    ),
    TestDataStructure(
        "alagly",
        {"-ipf", "fda.pfa", "-s", "conf.gro", "-frame", "0"},
        "",
        "result.pdb",
        "ref.pdb"
    ),
    TestDataStructure(
        "alagly",
        {"-ipf", "fda.pfa", "-s", "conf.gro", "-frame", "all", "-t", "1000"},
        "",
        "result2.pdb",
        "ref2.pdb"
    )
));
