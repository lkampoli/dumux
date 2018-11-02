dune_symlink_to_source_files(FILES params.input)

# isothermal tests
dune_add_test(NAME test_2pnc_fuelcell_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FuelCellBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pnc_fuelcell_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_fuelcell_box-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_fuelcell_box -ParameterFile params.input -Problem.Name test_2pnc_fuelcell_box")

dune_add_test(NAME test_2pnc_fuelcell_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FuelCellCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pnc_fuelcell_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_fuelcell_tpfa-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pnc_fuelcell_tpfa -ParameterFile params.input -Problem.Name test_2pnc_fuelcell_tpfa")


dune_add_test(NAME test_2pncni_fuelcell_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FuelCellNIBoxTypeTag NONISOTHERMAL=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pncni_fuelcell_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pncni_fuelcell_box-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncni_fuelcell_box -ParameterFile params.input -Problem.Name test_2pncni_fuelcell_box")

#install sources
install(FILES
main.cc
problem.hh
spatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pnc/fuelcell)
set(CMAKE_BUILD_TYPE Release)