if (NOT (TETGEN_INCLUDES AND TETGEN_LIBRARIES))
    find_path(TETGEN_INCLUDES NAMES tetgen.h HINTS ${TETGEN_INCLUDE_DIR} )
    find_library(TETGEN_LIBRARIES tet HINTS ${TETGEN_LIBRARY_DIR} )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    TETGEN DEFAULT_MSG TETGEN_INCLUDES TETGEN_LIBRARIES)

mark_as_advanced(TETGEN_INCLUDES TETGEN_LIBRARIES)
