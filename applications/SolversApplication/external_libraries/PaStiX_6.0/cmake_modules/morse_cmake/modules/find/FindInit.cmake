###
#
# @copyright (c) 2018 Inria. All rights reserved.
#
###
#
#  @file FindInit.cmake
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
#  @author Florent Pruvost
#  @date 24-04-2018
#
###


# This include is required to check symbols of libs
include(CheckFunctionExists)

# This include is required to check defines in headers
include(CheckIncludeFiles)

# Factorize some piece of code
include(FindCommon)

# To find headers and libs
include(FindHeadersAndLibs)

# To transform relative path into absolute for a list of libraries
include(LibrariesAbsolutePath)
include(FindPkgconfigLibrariesAbsolutePath)

# Some macros to print status when search for headers and libs
include(PrintFindStatus)

##
## @end file FindInit.cmake
##
