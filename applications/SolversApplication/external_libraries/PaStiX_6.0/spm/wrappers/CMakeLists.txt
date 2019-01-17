###
#
#  @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.0
#  @author Mathieu Faverge
#  @date 2017-05-22
#
###

if (SPM_WITH_FORTRAN)
  add_subdirectory( fortran90 )
endif()

if (BUILD_SHARED_LIBS)
  add_subdirectory( python )
else()
  message(STATUS "--- Python wrapper is disabled with static libraries")
endif()
