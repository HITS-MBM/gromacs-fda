set(testname "FDAGraphTest")
set(exename "fda-graph-test")

gmx_build_unit_test(
    ${testname}
    ${exename}
    # files with code for test fixtures
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
