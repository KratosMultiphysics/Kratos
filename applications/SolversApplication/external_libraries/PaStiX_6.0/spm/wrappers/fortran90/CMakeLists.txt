###
#
#  @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.0
#  @author Mathieu Faverge
#  @date 2017-05-22
#
###
cmake_minimum_required (VERSION 3.1)

# Coherce CMake to install the generated .mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

add_library( spmf
  src/spm_enums.F90
  src/spmf.f90 )

if ( SPM_INT64 )
  set_source_files_properties(
    src/spm_enums.F90
    PROPERTIES COMPILE_DEFINITIONS "SPM_INT_KIND=c_int64_t")
else()
  set_source_files_properties(
    src/spm_enums.F90
    PROPERTIES COMPILE_DEFINITIONS "SPM_INT_KIND=c_int32_t")
endif()

target_link_libraries( spmf spm )
install(TARGETS spmf
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib )

#
# Add examples
#
set (EXAMPLES
  spm_driver.f90
  spm_user.f90
  )

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} examples/${_file})
  target_link_libraries(${_name_we} spmf)

  install(TARGETS ${_name_we}       RUNTIME DESTINATION examples )
  install(FILES   examples/${_file}         DESTINATION examples )

  add_test(fortran_${_name_we} ./${_name_we})

endforeach()

