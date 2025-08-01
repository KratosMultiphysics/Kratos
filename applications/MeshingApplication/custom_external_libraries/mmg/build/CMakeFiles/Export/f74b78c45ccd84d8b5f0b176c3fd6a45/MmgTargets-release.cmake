#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Mmg::libmmg2d_so" for configuration "Release"
set_property(TARGET Mmg::libmmg2d_so APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Mmg::libmmg2d_so PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmmg2d.so.5.7.0"
  IMPORTED_SONAME_RELEASE "libmmg2d.so.5"
  )

list(APPEND _cmake_import_check_targets Mmg::libmmg2d_so )
list(APPEND _cmake_import_check_files_for_Mmg::libmmg2d_so "${_IMPORT_PREFIX}/lib/libmmg2d.so.5.7.0" )

# Import target "Mmg::libmmgs_so" for configuration "Release"
set_property(TARGET Mmg::libmmgs_so APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Mmg::libmmgs_so PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmmgs.so.5.7.0"
  IMPORTED_SONAME_RELEASE "libmmgs.so.5"
  )

list(APPEND _cmake_import_check_targets Mmg::libmmgs_so )
list(APPEND _cmake_import_check_files_for_Mmg::libmmgs_so "${_IMPORT_PREFIX}/lib/libmmgs.so.5.7.0" )

# Import target "Mmg::libmmg3d_so" for configuration "Release"
set_property(TARGET Mmg::libmmg3d_so APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Mmg::libmmg3d_so PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmmg3d.so.5.7.0"
  IMPORTED_SONAME_RELEASE "libmmg3d.so.5"
  )

list(APPEND _cmake_import_check_targets Mmg::libmmg3d_so )
list(APPEND _cmake_import_check_files_for_Mmg::libmmg3d_so "${_IMPORT_PREFIX}/lib/libmmg3d.so.5.7.0" )

# Import target "Mmg::libmmg_so" for configuration "Release"
set_property(TARGET Mmg::libmmg_so APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Mmg::libmmg_so PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmmg.so.5.7.0"
  IMPORTED_SONAME_RELEASE "libmmg.so.5"
  )

list(APPEND _cmake_import_check_targets Mmg::libmmg_so )
list(APPEND _cmake_import_check_files_for_Mmg::libmmg_so "${_IMPORT_PREFIX}/lib/libmmg.so.5.7.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
