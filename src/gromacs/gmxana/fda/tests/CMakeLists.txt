set(testname "FDAGraphTest")
set(exename "fda-graph-test")

gmx_add_gtest_executable(
    ${exename}
    FDABoostGraphTest.cpp
    FDAGetStressTest.cpp
    FDAGraphTest.cpp
    FDAShortestPathTest.cpp
    FDAViewStressTest.cpp
    PDBTest.cpp
)

gmx_register_integration_test(
    ${testname}
    ${exename}
)
