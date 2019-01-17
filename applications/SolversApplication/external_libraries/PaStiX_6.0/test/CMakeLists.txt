###
#
#  @copyright 2013-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.1
#  @author Mathieu Faverge
#  @date 2018-07-16
#
###
include(RulesPrecisions)

## Generate the test library files for all required precisions
set(LIB_SOURCES
  z_bcsc_tests.c
)

## reset variables
set(generated_libfiles "")

precisions_rules_py(generated_libfiles
   "${LIB_SOURCES}"
   PRECISIONS "p;s;d;c;z")

add_library(bcsc_test
   ${generated_libfiles}
   )

target_link_libraries( bcsc_test
  pastix
  pastix_kernels
  ${SPM_LIBRARY}
  ${LAPACKE_LIBRARIES_DEP}
  )

## Generate all test executables
set (TESTS
  bcsc_norm_tests.c
  bcsc_spmv_tests.c
  )

foreach (_file ${TESTS})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} pastix bcsc_test)
endforeach()

## CTest execution
set( BCSC_TESTS
  bcsc_norm_tests bcsc_spmv_tests )

# Laplacian
foreach(example ${BCSC_TESTS} )
  foreach(arithm ${PASTIX_PRECISIONS} )
    add_test(test_lap_${arithm}_${example} ./${example} --lap ${arithm}:10:10:10)
  endforeach()
endforeach()

# RSA
foreach(example ${BCSC_TESTS} )
  add_test(test_rsa_${example} ./${example} --rsa ${CMAKE_SOURCE_DIR}/test/matrix/small.rsa)
endforeach()

# Matrix Market
foreach(example ${BCSC_TESTS} )
  add_test(test_mm_${example} ./${example} --mm ${CMAKE_SOURCE_DIR}/test/matrix/young4c.mtx)
endforeach()

# Harwell Boeing
foreach(example ${BCSC_TESTS} )
  add_test(test_hb_${example} ./${example} --hb ${CMAKE_SOURCE_DIR}/test/matrix/orsirr.rua)
endforeach()

# Matrix Market - Hermitian
foreach(example ${BCSC_TESTS} )
  add_test(test_mm2_${example} ./${example} --mm ${CMAKE_SOURCE_DIR}/test/matrix/mhd1280b.mtx)
endforeach()

###################################
#   Low-rank testings
###

## Tests requires TMG lib to generate random matrices.
find_package(TMG)

if(TMG_FOUND)
  ## Generate the test library files for all required precisions
  set(LIB_SOURCES
    z_lowrank_tests.c
    z_bvec_tests.c
    )

  ## reset variables
  set(generated_libfiles "")

  precisions_rules_py(generated_libfiles
    "${LIB_SOURCES}"
    PRECISIONS "p;s;d;c;z")

  add_library(pastix_tests
    ${generated_libfiles}
    )

  target_link_libraries( pastix_tests
    pastix
    pastix_kernels
    ${LAPACKE_LIBRARIES_DEP}
    ${TMG_LIBRARIES_DEP}
    ${EXTRA_LIBS}
    m)

  set(TESTS_SOURCES
    z_ge2lr_tests.c
    z_rradd_tests.c
    z_lrmm_tests.c
    )

  ## reset variables
  set(generated_exfiles)

  precisions_rules_py(generated_exfiles
    "${TESTS_SOURCES}"
    PRECISIONS "p;s;d;c;z")

  set( testings
    ${generated_exfiles}
    bcsc_spmv_time.c
    bvec_gemv_tests.c
    bvec_tests.c )

  foreach (_file ${testings})
    get_filename_component(_name_we ${_file} NAME_WE)
    add_executable(${_name_we} ${_file})
    target_link_libraries(${_name_we}
      pastix_tests
      pastix_kernels
      pastix
      ${LAPACKE_LIBRARIES_DEP}
      ${TMG_LIBRARIES_DEP}
      ${EXTRA_LIBS}
      m)

    add_test(test_${_name_we} ./${_name_we} -9 1000:10:1 )
  endforeach()

endif()

