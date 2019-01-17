###
#
#  @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.1
#  @author Mathieu Faverge
#  @date 2018-07-16
#
###
cmake_minimum_required (VERSION 3.1)

# Coherce CMake to install the generated .mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

add_library( pastixf
  src/pastix_enums.F90
  src/pastixf.f90 )

if ( PASTIX_INT64 )
  set_source_files_properties(
    src/pastix_enums.F90
    PROPERTIES COMPILE_DEFINITIONS "PASTIX_INT_KIND=c_int64_t")
else()
  set_source_files_properties(
    src/pastix_enums.F90
    PROPERTIES COMPILE_DEFINITIONS "PASTIX_INT_KIND=c_int32_t")
endif()

target_link_libraries( pastixf pastix )
install(TARGETS pastixf
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib )

#
# Add examples
#
set (EXAMPLES
  fsimple.f90
  flaplacian.f90
  fstep-by-step.f90
  )

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} examples/${_file})
  target_link_libraries(${_name_we} pastixf spmf)

  install(TARGETS ${_name_we}       RUNTIME DESTINATION examples )
  install(FILES   examples/${_file}         DESTINATION examples )

  add_test(fortran_${_name_we} ./${_name_we})

endforeach()

#
# Add the fmultilap example using OpenMP
#
set( _file fmultilap.f90 )

get_filename_component(_name_we ${_file} NAME_WE)
add_executable(${_name_we} examples/${_file})
target_link_libraries(${_name_we} pastixf spmf)

install(TARGETS ${_name_we}       RUNTIME DESTINATION examples )
install(FILES   examples/${_file}         DESTINATION examples )
install(FILES   examples/test_seq.in      DESTINATION examples )
install(FILES   examples/test_mt.in       DESTINATION examples )

# Add OpenMP if available
find_package(OpenMP)
if (OpenMP_Fortran_FOUND)
  set_target_properties( fmultilap
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS} )
  target_link_libraries( fmultilap
    ${OpenMP_Fortran_LIBRARIES})
endif()

# Add specific test
include(FindUnixCommands)

set( _seq "${CMAKE_CURRENT_SOURCE_DIR}/examples/test_seq.in")
set( _mt  "${CMAKE_CURRENT_SOURCE_DIR}/examples/test_mt.in" )
set( _exe "${CMAKE_CURRENT_BINARY_DIR}/fmultilap"           )

if (BASH)
  add_test( NAME fortran_${_name_we}_seq COMMAND ${BASH} -c "${_exe} < ${_seq}" )
  add_test( NAME fortran_${_name_we}_mt  COMMAND ${BASH} -c "${_exe} < ${_mt}"  )
endif()

#
# End of fmultilap
#
