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

SET ( GET_MMG_TESTS "FALSE" )
SET ( GET_MMG2D_TESTS "FALSE" )
SET ( GET_MMGS_TESTS "FALSE" )
SET ( GET_MMG3D_TESTS "FALSE" )

# Check if the ci_tests directory exists
IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests )

  # First download of the tests
  MESSAGE("-- Create the continuous integration directory (ci_tests)"
    " and download the test cases. May be very long...")
  FILE(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests)

  SET ( GET_MMG_TESTS "TRUE" )
  SET ( GET_MMG2D_TESTS "TRUE" )
  SET ( GET_MMGS_TESTS "TRUE" )
  SET ( GET_MMG3D_TESTS "TRUE" )

ELSE ( )

  # Check if the tests are up to date
  #--------------> mmg
  IF ( EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg.version )
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmg.version OLD_MMG_MD5)
  ELSE ( )
    SET ( OLD_MMG_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmRktsVkFDTGlfdzQ
    ${CMAKE_SOURCE_DIR}/ci_tests/mmg.version
    STATUS MMG_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMG_VERSION_STATUS 0 MMG_VERSION_STATUS_0)
  LIST(GET MMG_VERSION_STATUS 1 MMG_VERSION_STATUS_1)

  IF ( MMG_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmg.version MMG_MD5)

    IF ( NOT (${OLD_MMG_MD5} MATCHES ${MMG_MD5}) )
      SET ( GET_MMG_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMG_VERSION_STATUS_1}.
 Try to get it at the following link:
        https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmRktsVkFDTGlfdzQ
 then untar it in the project directory (mmg/ by default).")
  ENDIF()

  #--------------> mmg2d
  IF ( EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg2d.version )
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmg2d.version OLD_MMG2D_MD5)
  ELSE ( )
    SET ( OLD_MMG2D_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmV3BlUER4M0Z4MGs
    ${CMAKE_SOURCE_DIR}/ci_tests/mmg2d.version
    STATUS MMG2D_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMG2D_VERSION_STATUS 0 MMG2D_VERSION_STATUS_0)
  LIST(GET MMG2D_VERSION_STATUS 1 MMG2D_VERSION_STATUS_1)

  IF ( MMG2D_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmg2d.version MMG2D_MD5)

    IF ( NOT (${OLD_MMG2D_MD5} MATCHES ${MMG2D_MD5}) )
      SET ( GET_MMG2D_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMG2D_VERSION_STATUS_1}.
 Try to get it at the following link:
        https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmV3BlUER4M0Z4MGs
 then untar it in the project directory (mmg/ by default).")
  ENDIF()

    #--------------> mmgs
  IF ( EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version )
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version OLD_MMGS_MD5)
  ELSE ( )
    SET ( OLD_MMGS_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmSWtGV295a28wU2c
    ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version
    STATUS MMGS_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMGS_VERSION_STATUS 0 MMGS_VERSION_STATUS_0)
  LIST(GET MMGS_VERSION_STATUS 1 MMGS_VERSION_STATUS_1)

  IF ( MMGS_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.version MMGS_MD5)

    IF ( NOT (${OLD_MMGS_MD5} MATCHES ${MMGS_MD5}) )
      SET ( GET_MMGS_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMGS_VERSION_STATUS_1}.
 Try to get it at the following link:
       https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmSWtGV295a28wU2c
 then untar it in the project directory (mmg/ by default).")
  ENDIF()

  #--------------> mmg3d
  IF ( EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version )
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version OLD_MMG3D_MD5)
  ELSE ( )
    SET ( OLD_MMG3D_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmSWtGV295a28wU2c
    ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version
    STATUS MMG3D_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMG3D_VERSION_STATUS 0 MMG3D_VERSION_STATUS_0)
  LIST(GET MMG3D_VERSION_STATUS 1 MMG3D_VERSION_STATUS_1)

  IF ( MMG3D_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.version MMG3D_MD5)

    IF (NOT (${OLD_MMG3D_MD5} MATCHES ${MMG3D_MD5}))
      SET ( GET_MMG3D_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMG3D_VERSION_STATUS_1}.
 Try to get it at the following link:
        https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmSWtGV295a28wU2c
 then untar it in the project directory (mmg/ by default).")
  ENDIF()

ENDIF()

# Download the tests if needed
#--------------> mmg
IF ( GET_MMG_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmg test cases. May be very long...")
  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmdXVkS3QzdWEtdFU
    ${CMAKE_SOURCE_DIR}/ci_tests/mmg.tgz
    SHOW_PROGRESS)
  IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg.tgz )
    MESSAGE("\n")
    MESSAGE(WARNING "Fail to automatically download the mmg test cases.
Try to get it at the following link:
       https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmdXVkS3QzdWEtdFU
then untar it in the project directory (mmg/ by default).")
  ELSE()
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar xzf
      ${CMAKE_SOURCE_DIR}/ci_tests/mmg.tgz
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests/
      )
    IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg.tgz )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmg "
        "test cases directory (mmg.tgz).
Try to untar it by hand in the project directory"
        " (mmg/ci_tests/ by default).")
    ENDIF()
    FILE(REMOVE ${CMAKE_SOURCE_DIR}/ci_tests/mmg.tgz)
  ENDIF ()
ENDIF ()


#--------------> mmg2d
IF ( GET_MMG2D_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmg2d test cases. May be very long...")
  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmX0hrRWJWTDBETHc
    ${CMAKE_SOURCE_DIR}/ci_tests/mmg2d.tgz
    SHOW_PROGRESS)
  IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg2d.tgz )
    MESSAGE("\n")
    MESSAGE(WARNING "Fail to automatically download the mmg2d test cases.
Try to get it at the following link:
       https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmX0hrRWJWTDBETHc
then untar it in the project directory (mmg/ by default).")
  ELSE()
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar xzf
      ${CMAKE_SOURCE_DIR}/ci_tests/mmg2d.tgz
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests/
      )
    IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg2d.tgz )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmg2d "
        "test cases directory (mmg2d.tgz).
Try to untar it by hand in the project directory"
        " (mmg/ci_tests/ by default).")
    ENDIF()
    FILE(REMOVE ${CMAKE_SOURCE_DIR}/ci_tests/mmg2d.tgz)
  ENDIF ()
ENDIF ()

#--------------> mmgs
IF ( GET_MMGS_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmgs test cases. May be very long...")
  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmcVdZb1EzaTR3ZlU
    ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz
    SHOW_PROGRESS)
  IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz )
    MESSAGE("\n")
    MESSAGE(WARNING "Fail to automatically download the mmgs test cases.
Try to get it at the following link:
       https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmcVdZb1EzaTR3ZlU
then untar it in the project directory (mmg/ by default).")
  ELSE()
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar xzf
      ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests/
      )
    IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmgs "
        "test cases directory (mmgs.tgz).
Try to untar it by hand in the project directory"
        " (mmg/ci_tests/ by default).")
    ENDIF()
    FILE(REMOVE ${CMAKE_SOURCE_DIR}/ci_tests/mmgs.tgz)
  ENDIF ()
ENDIF ()

#--------------> mmg3d
SET(ADDRESS
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmWGxhMnAzMGFrNTg
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmclVnbTRqUXVfNUE
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmdWNLTFh1THl1U0k
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmX1d2WGJkaHUxV1E
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmbmY1R0EyelhETW8
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmd3lnTTZQRTZxLW8
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmRHk0ZTREdzJpYXc
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmd09YVmtEcXZIWFU
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmLTFXUGhGTVY2dEE
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmUGI3UU1UMExJTTQ
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmSWZJTEt1aGVCZDA
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmLUZqNjhwajBKZWM
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmZG9lUTdGa1d5ZjA
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmU2l0N0N4X0lyRkE
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmcmVnY2NsSzBxVkk
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmRWIxRS1yMUJyZkk
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmM05MVEhvZzdHckE
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmbFhOaE5UVlNQOFk
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmLWp4SFpyN3c4d0U
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmZm0tODJtRXM2eHc
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmZ29NVVpyRlFob1k
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmd0lmd1pHR1VKUlU
  https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmY0VWR2pCbXRWLWM
  )

SET(FILENAME
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.aa
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ab
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ac
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ad
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ae
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.af
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ag
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ah
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ai
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.aj
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ak
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.al
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.am
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.an
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ao
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ap
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.aq
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.ar
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.as
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.at
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.au
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.av
  ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d.tgz.aw
  )


IF ( GET_MMG3D_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmg3d test cases. May be very long...")

  SET(LOAD_OK 1)

  FOREACH( i RANGE 0 22)
    LIST(GET ADDRESS  ${i} ADDRESS_i)
    LIST(GET FILENAME ${i} FILENAME_i)

    FILE(DOWNLOAD ${ADDRESS_i}
      ${FILENAME_i}
      SHOW_PROGRESS)
    IF ( NOT EXISTS ${FILENAME_i} )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically download the mmg3d test cases
Try to get it at the following link:
       ${ADDRESS_i}
then untar it in the project directory mmg/ by default).")
      SET ( LOAD_OK 0 )
      BREAK()
    ENDIF()

  ENDFOREACH()

  IF ( ${LOAD_OK} )
    EXECUTE_PROCESS(
      COMMAND cat ${FILENAME}
      COMMAND tar -xzf -
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ci_tests/
      TIMEOUT 10000
      )
    IF ( NOT EXISTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmg3d"
        "test cases directory (mmg3d.tgz.*).
Try to untar it by hand in the project directory"
        " (mmg/ci_tests/ by default):
        cat mmg3d.tgz.* | tar xzvf - ")
    ENDIF()

    FILE(REMOVE ${FILENAME})
  ENDIF()

ENDIF()
