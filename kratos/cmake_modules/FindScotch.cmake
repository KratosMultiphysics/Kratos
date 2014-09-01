# Pastix requires SCOTCH or METIS (partitioning and reordering tools)

if (NOT (SCOTCH_INCLUDES AND SCOTCH_LIBRARIES))
    find_path(SCOTCH_INCLUDES
        NAMES scotch.h
        PATHS $ENV{SCOTCHDIR}
        PATH_SUFFIXES scotch
        )


    find_library(SCOTCH_LIBRARY      scotch      PATHS $ENV{SCOTCHDIR})
    find_library(PTSCOTCH_LIBRARY    ptscotch    PATHS $ENV{SCOTCHDIR})
    find_library(SCOTCHERR_LIBRARY   scotcherr   PATHS $ENV{SCOTCHDIR})
    find_library(PTSCOTCHERR_LIBRARY ptscotcherr PATHS $ENV{SCOTCHDIR})

    set(SCOTCH_LIBRARIES
        "${SCOTCH_LIBRARY};${PTSCOTCH_LIBRARY};${SCOTCHERR_LIBRARY};${PTSCOTCHERR_LIBRARY}"
        )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    SCOTCH DEFAULT_MSG SCOTCH_INCLUDES SCOTCH_LIBRARIES)

mark_as_advanced(SCOTCH_INCLUDES SCOTCH_LIBRARIES)
