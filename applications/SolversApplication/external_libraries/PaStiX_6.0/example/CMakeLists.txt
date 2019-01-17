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
include(GenPkgConfig)

add_subdirectory(old)

set (INTERNAL_EXAMPLES
  dump_rank.c
  ordering_grid.c
  )

foreach (_file ${INTERNAL_EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} pastix)
  target_link_libraries(${_name_we} ${BLAS_LIBRARIES})
endforeach()

set (EXAMPLES
  analyze.c
  bench_facto.c
  compress.c
  multidof.c
  personal.c
  reentrant.c
  refinement.c
  schur.c
  simple.c
  step-by-step.c
  )

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} pastix)
  target_link_libraries(${_name_we} ${BLAS_LIBRARIES})

  install(TARGETS ${_name_we} RUNTIME DESTINATION examples )
  install(FILES   ${_file}            DESTINATION examples )
endforeach()

### Install a makefile for user compilation test
set( EXAMPLE_LIBS ${DEPS_LIBRARIES} )
if ( EXAMPLE_LIBS )
  gpc_convert_libstyle_to_pkgconfig( EXAMPLE_LIBS )
  string(REPLACE ";" " " EXAMPLE_LIBS "${EXAMPLE_LIBS}")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Makefile.in"
  "${CMAKE_CURRENT_BINARY_DIR}/make/Makefile" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/make/Makefile"
  DESTINATION examples)

include(CTestLists.txt)
