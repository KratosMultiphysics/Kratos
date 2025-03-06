#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "OpenSubdiv::osdCPU_static" for configuration ""
set_property(TARGET OpenSubdiv::osdCPU_static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(OpenSubdiv::osdCPU_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libosdCPU.a"
  )

list(APPEND _cmake_import_check_targets OpenSubdiv::osdCPU_static )
list(APPEND _cmake_import_check_files_for_OpenSubdiv::osdCPU_static "${_IMPORT_PREFIX}/lib/libosdCPU.a" )

# Import target "OpenSubdiv::osdGPU_static" for configuration ""
set_property(TARGET OpenSubdiv::osdGPU_static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(OpenSubdiv::osdGPU_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libosdGPU.a"
  )

list(APPEND _cmake_import_check_targets OpenSubdiv::osdGPU_static )
list(APPEND _cmake_import_check_files_for_OpenSubdiv::osdGPU_static "${_IMPORT_PREFIX}/lib/libosdGPU.a" )

# Import target "OpenSubdiv::osdCPU" for configuration ""
set_property(TARGET OpenSubdiv::osdCPU APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(OpenSubdiv::osdCPU PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libosdCPU.so.3.6.0"
  IMPORTED_SONAME_NOCONFIG "libosdCPU.so.3.6.0"
  )

list(APPEND _cmake_import_check_targets OpenSubdiv::osdCPU )
list(APPEND _cmake_import_check_files_for_OpenSubdiv::osdCPU "${_IMPORT_PREFIX}/lib/libosdCPU.so.3.6.0" )

# Import target "OpenSubdiv::osdGPU" for configuration ""
set_property(TARGET OpenSubdiv::osdGPU APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(OpenSubdiv::osdGPU PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libosdGPU.so.3.6.0"
  IMPORTED_SONAME_NOCONFIG "libosdGPU.so.3.6.0"
  )

list(APPEND _cmake_import_check_targets OpenSubdiv::osdGPU )
list(APPEND _cmake_import_check_files_for_OpenSubdiv::osdGPU "${_IMPORT_PREFIX}/lib/libosdGPU.so.3.6.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
