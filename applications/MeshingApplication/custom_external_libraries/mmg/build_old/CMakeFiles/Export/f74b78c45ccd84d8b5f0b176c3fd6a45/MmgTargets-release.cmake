#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Mmg::libmmg2d_a" for configuration "Release"
set_property(TARGET Mmg::libmmg2d_a APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Mmg::libmmg2d_a PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmmg2d.a"
  )

list(APPEND _cmake_import_check_targets Mmg::libmmg2d_a )
list(APPEND _cmake_import_check_files_for_Mmg::libmmg2d_a "${_IMPORT_PREFIX}/lib/libmmg2d.a" )

# Import target "Mmg::libmmgs_a" for configuration "Release"
set_property(TARGET Mmg::libmmgs_a APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Mmg::libmmgs_a PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmmgs.a"
  )

list(APPEND _cmake_import_check_targets Mmg::libmmgs_a )
list(APPEND _cmake_import_check_files_for_Mmg::libmmgs_a "${_IMPORT_PREFIX}/lib/libmmgs.a" )

# Import target "Mmg::libmmg3d_a" for configuration "Release"
set_property(TARGET Mmg::libmmg3d_a APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Mmg::libmmg3d_a PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmmg3d.a"
  )

list(APPEND _cmake_import_check_targets Mmg::libmmg3d_a )
list(APPEND _cmake_import_check_files_for_Mmg::libmmg3d_a "${_IMPORT_PREFIX}/lib/libmmg3d.a" )

# Import target "Mmg::libmmg_a" for configuration "Release"
set_property(TARGET Mmg::libmmg_a APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Mmg::libmmg_a PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmmg.a"
  )

list(APPEND _cmake_import_check_targets Mmg::libmmg_a )
list(APPEND _cmake_import_check_files_for_Mmg::libmmg_a "${_IMPORT_PREFIX}/lib/libmmg.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
