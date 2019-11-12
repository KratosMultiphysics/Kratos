#          __        ____ ____
# ___  ___/  |_ __ _/_   /_   |
# \  \/ /\   __\  |  \   ||   |
#  \   /  |  | |  |  /   ||   |
#   \_/   |__| |____/|___||___|
#
#  License: BSD License ; see LICENSE
#


cmake_minimum_required( VERSION 3.4.3 )

project( vtu11 CXX )

set( CMAKE_CXX_STANDARD 11 )

# Directories to include header files from
include_directories( vtu11 )
include_directories( ${CMAKE_CURRENT_BINARY_DIR}/config )

# -------------------- Set up testrunner --------------------------

OPTION ( VTU11_ENABLE_TESTS "Enable unit tests." ON )

if( VTU11_ENABLE_TESTS )

    file( GLOB TEST_FILES test/*.cpp)

    # Add an executable and link to the library
    add_executable( vtu11_testrunner ${TEST_FILES} )

    file(COPY test/testfiles DESTINATION . )

endif( VTU11_ENABLE_TESTS )

# -------------------- propagate include list --------------------------

get_directory_property( hasParent PARENT_DIRECTORY )

if(hasParent)

    # Write all include directories into ALL_INCLUDES variable
    get_property( ALL_INCLUDES DIRECTORY PROPERTY INCLUDE_DIRECTORIES )

    # Set new variable VTU11_INCLUDE_DIRECTORIES also in parent scope
    set( VTU11_INCLUDE_DIRECTORIES ${ALL_INCLUDES} PARENT_SCOPE )

endif()

# ---------------------------------------------------------------------

OPTION ( VTU11_STRICT_COMPILER "Enable more warnings as errors." OFF )

if( VTU11_STRICT_COMPILER )

    # If compiler is g++: Enable further  warnings and treat all warnings as errors
    if( CMAKE_COMPILER_IS_GNUCXX )

        set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -pedantic -Wall -Wextra -Wcast-align -Wsuggest-attribute=cold       \
            -Wsuggest-attribute=pure -Wimport -Wsuggest-final-methods -Wsuggest-attribute=const -Wsuggest-attribute=format \
            -Wsuggest-attribute=malloc -Wsuggest-attribute=noreturn -Wformat-y2k -Wpacked -Wno-restrict -Wpadded           \
            -Wswitch-enum -Wwrite-strings -Wformat-nonliteral -Wformat-security -Wcast-qual -Wsuggest-override    \
            -Wsuggest-final-types -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs       \
            -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo     \
            -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror" )

    endif( CMAKE_COMPILER_IS_GNUCXX )

endif( VTU11_STRICT_COMPILER )

# -------------------- zlib dependency ---------------------------

find_package( ZLIB )

if ( ZLIB_FOUND )

    message( STATUS "Enabling vtu11 with zlib compression" )

    include_directories( ${ZLIB_INCLUDE_DIRS} )

    if( VTU11_ENABLE_TESTS )
        target_link_libraries( vtu11_testrunner ${ZLIB_LIBRARIES} )
    endif( VTU11_ENABLE_TESTS )
    
else ( ZLIB_FOUND )

    message( STATUS "Enabling vtu11 without zlib compression" )
    
endif( ZLIB_FOUND )

# ---------------------- Copy / install ----------------------------

CONFIGURE_FILE ( vtu11_config.hpp.in config/inc/vtu11_config.hpp )
