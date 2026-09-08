# FindTBB.cmake
# Finds the Intel Threading Building Blocks (TBB) library and creates the
# imported target TBB::tbb used by Kratos when KRATOS_USE_TBB is ON.
#
# This module sets the following variables:
#   TBB_FOUND        - TRUE if TBB headers and library were found
#   TBB_INCLUDE_DIRS - Directory containing TBB headers
#   TBB_LIBRARIES    - Full path to the TBB library
#
# And defines the imported target:
#   TBB::tbb         - The TBB library target

# Search for the TBB header
find_path(TBB_INCLUDE_DIR
    NAMES tbb/tbb.h
    HINTS
        ${TBB_DIR}
        ${TBB_ROOT}
        $ENV{TBB_DIR}
        $ENV{TBB_ROOT}
        $ENV{TBBROOT}
    PATHS
        /usr/include
        /usr/local/include
        /opt/intel/tbb/include
        /opt/tbb/include
)

# Search for the TBB library
find_library(TBB_LIBRARY
    NAMES tbb
    HINTS
        ${TBB_DIR}
        ${TBB_ROOT}
        $ENV{TBB_DIR}
        $ENV{TBB_ROOT}
        $ENV{TBBROOT}
    PATH_SUFFIXES
        lib
        lib64
        lib/intel64
        lib/intel64/gcc4.8
    PATHS
        /usr/lib
        /usr/lib64
        /usr/local/lib
        /usr/local/lib64
        /opt/intel/tbb/lib
        /opt/tbb/lib
)

find_library(TBB_MALLOC_LIBRARY
    NAMES tbbmalloc
    HINTS
        ${TBB_DIR}
        ${TBB_ROOT}
        $ENV{TBB_DIR}
        $ENV{TBB_ROOT}
        $ENV{TBBROOT}
    PATH_SUFFIXES
        lib
        lib64
        lib/intel64
        lib/intel64/gcc4.8
    PATHS
        /usr/lib
        /usr/lib64
        /usr/local/lib
        /usr/local/lib64
        /opt/intel/tbb/lib
        /opt/tbb/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TBB
    REQUIRED_VARS TBB_LIBRARY TBB_INCLUDE_DIR
)

if(TBB_FOUND)
    set(TBB_INCLUDE_DIRS ${TBB_INCLUDE_DIR})
    set(TBB_LIBRARIES    ${TBB_LIBRARY})

    if(NOT TARGET TBB::tbb)
        add_library(TBB::tbb UNKNOWN IMPORTED)
        set_target_properties(TBB::tbb PROPERTIES
            IMPORTED_LOCATION             "${TBB_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${TBB_INCLUDE_DIR}"
            INTERFACE_COMPILE_DEFINITIONS "KRATOS_USE_TBB"
        )
    endif()

    if(TBB_MALLOC_LIBRARY AND NOT TARGET TBB::tbbmalloc)
        add_library(TBB::tbbmalloc UNKNOWN IMPORTED)
        set_target_properties(TBB::tbbmalloc PROPERTIES
            IMPORTED_LOCATION             "${TBB_MALLOC_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${TBB_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(TBB_INCLUDE_DIR TBB_LIBRARY TBB_MALLOC_LIBRARY)
