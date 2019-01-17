###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2018 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
#  @file MorseInit.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 1.0.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 13-07-2012
#
###

# Path to Morse modules
get_filename_component(MORSE_CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_FILE} DIRECTORY CACHE)

# Global Morse options
option(MORSE_ENABLE_WARNING       "Enable warning messages" OFF)
option(MORSE_ENABLE_COVERAGE      "Enable flags for coverage test" OFF)
option(MORSE_ENABLE_COLOR_MESSAGE "Enable colors in messages" OFF)
#option(MORSE_VERBOSE_FIND_PACKAGE "Add additional messages concerning packages not found" OFF)
#message(STATUS "MORSE_VERBOSE_FIND_PACKAGE is set to OFF, turn it ON to get"
#        "   information about packages not found")


# This include is required to check symbols of libs in the main CMakeLists.txt
include(CheckFunctionExists)

# This include is required to check defines in headers
include(CheckIncludeFiles)

if (MORSE_ENABLE_COLOR_MESSAGE)
  # colorize messages
  include(ColorizeMessage)
endif()

# Define some auxilary flags
include(AuxilaryFlags)

# Define some variables to et info about ressources
include(Ressources)

# Add the path where we handle our FindFOO.cmake to seek for liraries
list(APPEND CMAKE_MODULE_PATH ${MORSE_CMAKE_MODULE_PATH}/find)

# To load some macros used in Finds (could be useful for other projects)
include(FindInit)

##
## @end file MorseInit.cmake
##
