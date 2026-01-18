## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
##
##  mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the mmg distribution only if you accept them.
## =============================================================================

###############################################################################
#####
#####         Download test cases for Mmg<x> code
#####
###############################################################################
OPTION ( CI_CONTEXT "Disable network test and download progress bar for CI runs" OFF)

# Continuous integration tests are hosted on the sync server of inria
# (sync.bordeaux.inria.fr) in the /web/users/html/mmg/ who is published in
# https://static.bordeaux.inria.fr/mmg/.
#
# To add new tests, you will have:
#  1. to be affiliated to Inria;
#  2. to ask for an account on the server (https://doc-si.inria.fr/display/SU/Sauvegardes#)
#  3. to be member of the sedbor group (or to ask to change the directory group)
#  3. to follow the /web/users/html/mmg/README advices.

MACRO ( DOWNLOAD_TESTS x )

  # Check if the tests are up to date
  SET ( OLD_MMG_MD5 "0" )
  IF ( EXISTS ${CI_DIR}/mmg${x}.version )
    FILE(MD5 ${CI_DIR}/mmg${x}.version OLD_MMG_MD5)
  ENDIF ( )

  FILE(DOWNLOAD https://static.bordeaux.inria.fr/mmg/mmg${x}.version
    ${CI_DIR}/mmg${x}.version
    STATUS MMG_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMG_VERSION_STATUS 0 MMG_VERSION_STATUS_0)
  LIST(GET MMG_VERSION_STATUS 1 MMG_VERSION_STATUS_1)

  SET ( GET_MMG_TESTS "FALSE" )
  IF ( MMG_VERSION_STATUS_0 MATCHES 0)

    FILE(MD5 ${CI_DIR}/mmg${x}.version MMG_MD5)
    IF ( NOT (${OLD_MMG_MD5} MATCHES ${MMG_MD5}) )
      SET ( GET_MMG_TESTS "TRUE" )
    ENDIF ()

  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMG_VERSION_STATUS_1}.
 Try to get the test cases at the following link:
        https://static.bordeaux.inria.fr/mmg/mmg${x}.tgz
 Then untar it in ${CI_DIR}/mmg${x}.")

  ENDIF()

  # Get tests
  IF ( GET_MMG_TESTS MATCHES "TRUE" )
    MESSAGE("-- Mmg${x} test cases download. May take a while...")

    IF ( NOT ${CI_CONTEXT} )
      FILE(DOWNLOAD https://static.bordeaux.inria.fr/mmg/mmg${x}.tgz
        ${CI_DIR}/mmg${x}.tgz
        SHOW_PROGRESS)
    ELSE()
      FILE(DOWNLOAD https://static.bordeaux.inria.fr/mmg/mmg${x}.tgz
        ${CI_DIR}/mmg${x}.tgz)
    ENDIF()

    IF ( NOT EXISTS ${CI_DIR}/mmg${x}.tgz )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically download the mmg${x} test cases.
Try to get it at the following link:
       https://static.bordeaux.inria.fr/mmg/mmg${x}.tgz
then untar it in the ${CI_DIR}/mmg${x} directory.")
    ELSE()
      EXECUTE_PROCESS(
        COMMAND ${CMAKE_COMMAND} -E tar xzf ${CI_DIR}/mmg${x}.tgz
        WORKING_DIRECTORY ${CI_DIR}
        #RESULT_VARIABLE toto
        #COMMAND_ECHO STDOUT
        )
      IF ( NOT EXISTS ${CI_DIR}/mmg${x} )
        MESSAGE("\n")
        MESSAGE(WARNING "Fail to automatically untar the mmg${x} "
          "test cases directory (mmg${x}.tgz).
Try to untar it by hand in the ${CI_DIR} directory.")
      ELSE()
        FILE(REMOVE ${CI_DIR}/mmg${x}.tgz)
      ENDIF ()

    ENDIF()
  ENDIF()

ENDMACRO ( )

###############################################################################
#####
#####         Download test cases depending on user options
#####
###############################################################################
SET ( NO_CONNECTION 0 )

IF ( NOT ${CI_CONTEXT} )
  if ( MSVC )
    execute_process(
      COMMAND ping www.mmgtools.org -n 2
      OUTPUT_QUIET
      ERROR_QUIET
      RESULT_VARIABLE NO_CONNECTION
      )
  else ( )
    execute_process(
      COMMAND ping www.mmgtools.org -c 2
      OUTPUT_QUIET
      ERROR_QUIET
      RESULT_VARIABLE NO_CONNECTION
      )
  endif ( )

ENDIF()
if ( NOT (NO_CONNECTION EQUAL 0 OR NO_CONNECTION EQUAL 2) )
    set ( CONNECTED OFF )
    message ( STATUS "Offline mode: requires already downloaded test cases")
else()
  set ( CONNECTED ON )
endif()

IF ( CONNECTED )
  # Check if the ci_tests directory exists
  IF ( NOT EXISTS ${CI_DIR} )

    # First download of the tests
    MESSAGE("-- Creation of continuous integration directory: ${CI_DIR}")
    FILE(MAKE_DIRECTORY ${CI_DIR})

  ENDIF ( )

  #
  # Download the tests if needed
  ##--------------> mmg
  DOWNLOAD_TESTS ( "2d" )
  DOWNLOAD_TESTS ( "3d" )
  DOWNLOAD_TESTS ( "s" )
  DOWNLOAD_TESTS ( "" )
ENDIF ( )
