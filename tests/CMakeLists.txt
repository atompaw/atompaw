# Regression tests, run with ctest.
#
# The examples are run and their energies compared with the summary files
# shipped in example/; see run_example.sh for what is compared and the
# tolerances. The Hartree-Fock examples are not included: their inputs use
# "loggrid", whose grid has changed since their reference outputs were made,
# and they need the fixes of issue #12. example/Au/dirac has no reference.

set(ATP_TEST_DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/run_example.sh)

# atompaw_example(NAME DIR INPUT SUMMARY): run example/DIR/INPUT and compare
# the summary file SUMMARY it writes with example/DIR/SUMMARY
function(atompaw_example name dir input summary)
  add_test(NAME example.${name}
           COMMAND sh ${ATP_TEST_DRIVER} $<TARGET_FILE:atompaw>
                   ${CMAKE_CURRENT_BINARY_DIR}/${name}
                   ${PROJECT_SOURCE_DIR}/example/${dir}/${input}
                   ${summary}
                   ${PROJECT_SOURCE_DIR}/example/${dir}/${summary})
  set_tests_properties(example.${name} PROPERTIES LABELS "example" TIMEOUT 600)
endfunction()

atompaw_example(Li_lda Li/lda        Li.input Li)
atompaw_example(F_lda  F/lda         F.input  F)
atompaw_example(Na_lda Na/outputdata in       Na)

# All-electron carbon against the NIST atomic reference data for LDA (VWN),
# E = -37.425749 Ha = -74.851498 Ry, given to 1e-6 Ha; VWN needs libxc
if (LIBXC_OK)
  add_test(NAME nist.C_lda_vwn
           COMMAND sh ${ATP_TEST_DRIVER} $<TARGET_FILE:atompaw_aeonly>
                   ${CMAKE_CURRENT_BINARY_DIR}/nist_C_lda_vwn
                   ${CMAKE_CURRENT_SOURCE_DIR}/C_lda_vwn.input
                   C.AEO
                   =-74.851498 2e-6)
  set_tests_properties(nist.C_lda_vwn PROPERTIES LABELS "nist" TIMEOUT 600)
endif()
