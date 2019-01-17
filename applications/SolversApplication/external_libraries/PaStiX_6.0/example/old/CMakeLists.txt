###
#
#  @copyright 2013-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.1
#  @author Mathias Hastaran
#  @date 2018-07-16
#
###

set (EXAMPLES
  simple.c
  step-by-step.c
  )

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable("old_${_name_we}" ${_file})
  target_link_libraries("old_${_name_we}" pastix)

  install(TARGETS old_${_name_we} RUNTIME DESTINATION examples/old )
  install(FILES   ${_file}                DESTINATION examples/old )
endforeach()

### Install a makefile for user compilation test
set( EXAMPLE_LIBS
  ${DEPS_LIBRARIES}
  )
list(REMOVE_DUPLICATES EXAMPLE_LIBS)
string(REPLACE ";" " " EXAMPLE_LIBS "${EXAMPLE_LIBS}")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Makefile.in"
  "${CMAKE_CURRENT_BINARY_DIR}/make/Makefile" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/make/Makefile"
  DESTINATION examples/old )

