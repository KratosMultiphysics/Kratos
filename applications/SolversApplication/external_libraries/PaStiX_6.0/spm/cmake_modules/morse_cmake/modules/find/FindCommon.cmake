###
#
# @copyright (c) 2018 Inria. All rights reserved.
#
###
#
#  @file FindCommon.cmake
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
#  @date 13-04-2018
#
###

# clean these variables before using them in CMAKE_REQUIRED_* variables in
# check_function_exists
macro(finds_remove_duplicates)
  if (REQUIRED_DEFINITIONS)
    list(REMOVE_DUPLICATES REQUIRED_DEFINITIONS)
  endif()
  if (REQUIRED_INCDIRS)
    list(REMOVE_DUPLICATES REQUIRED_INCDIRS)
  endif()
  if (REQUIRED_FLAGS)
    list(REMOVE_DUPLICATES REQUIRED_FLAGS)
  endif()
  if (REQUIRED_LDFLAGS)
    list(REMOVE_DUPLICATES REQUIRED_LDFLAGS)
  endif()
  if (REQUIRED_LIBS)
    list(REVERSE REQUIRED_LIBS)
    list(REMOVE_DUPLICATES REQUIRED_LIBS)
    list(REVERSE REQUIRED_LIBS)
  endif()
endmacro()

##
## @end file FindCommon
##
