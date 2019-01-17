###
#
#  @copyright 2013-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @date 2013-06-24
#
###
include(RulesPrecisions)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

## Generate the test library files for all required precisions
set(LIB_SOURCES
  z_spm_tests.c
  core_zgeadd.c
  core_zplrnt.c
)

## reset variables
set(generated_libfiles "")

precisions_rules_py(generated_libfiles
   "${LIB_SOURCES}"
   PRECISIONS "p;s;d;c;z")

add_library(spm_test
   ${generated_libfiles}
   p_spm_tests.c
   get_options.c
   )

target_link_libraries( spm_test
  spm
  ${LAPACKE_LIBRARIES_DEP}
  ${CBLAS_LIBRARIES_DEP} )

## Generate all test executables
set (TESTS
  spm_convert_tests.c
  spm_norm_tests.c
  spm_matvec_tests.c
  spm_dof_expand_tests.c
  spm_dof_norm_tests.c
  spm_dof_matvec_tests.c
  )

foreach (_file ${TESTS})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} spm spm_test)
endforeach()

## CTest execution
set( SPM_TESTS
  spm_convert_tests spm_norm_tests spm_matvec_tests )
set( SPM_DOF_TESTS
  spm_dof_expand_tests spm_dof_norm_tests spm_dof_matvec_tests)

# Laplacian
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} )
  foreach(arithm ${RP_SPM_PRECISIONS} )
    add_test(test_lap_${arithm}_${example} ./${example}  --lap  ${arithm}:10:10:10:10.:2.)
    add_test(test_xlap_${arithm}_${example} ./${example} --xlap ${arithm}:6:10:12:5.:0.33)
  endforeach()
endforeach()

# RSA
foreach(example ${SPM_TESTS} )
  add_test(test_rsa_${example} ./${example} --rsa ${CMAKE_SOURCE_DIR}/tests/matrix/small.rsa)
endforeach()

# Matrix Market
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} )
  add_test(test_mm_${example} ./${example} --mm ${CMAKE_SOURCE_DIR}/tests/matrix/young4c.mtx)
endforeach()

# Harwell Boeing
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} )
  add_test(test_hb_${example} ./${example} --hb ${CMAKE_SOURCE_DIR}/tests/matrix/orsirr.rua)
endforeach()

# Matrix Market - Hermitian
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} )
  add_test(test_mm2_${example} ./${example} --mm ${CMAKE_SOURCE_DIR}/tests/matrix/mhd1280b.mtx)
endforeach()
