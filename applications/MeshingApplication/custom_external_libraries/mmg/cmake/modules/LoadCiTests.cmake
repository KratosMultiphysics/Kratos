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
IF ( NOT EXISTS ${CI_DIR} )

  # First download of the tests
  MESSAGE("-- Create the continuous integration directory ${CI_DIR}"
    " and download the test cases. May be very long...")
  FILE(MAKE_DIRECTORY ${CI_DIR})

  SET ( GET_MMG_TESTS "TRUE" )
  SET ( GET_MMG2D_TESTS "TRUE" )
  SET ( GET_MMGS_TESTS "TRUE" )
  SET ( GET_MMG3D_TESTS "TRUE" )

ELSE ( )

  # Check if the tests are up to date
  #--------------> mmg
  IF ( EXISTS ${CI_DIR}/mmg.version )
    FILE(MD5 ${CI_DIR}/mmg.version OLD_MMG_MD5)
  ELSE ( )
    SET ( OLD_MMG_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmRktsVkFDTGlfdzQ
    ${CI_DIR}/mmg.version
    STATUS MMG_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMG_VERSION_STATUS 0 MMG_VERSION_STATUS_0)
  LIST(GET MMG_VERSION_STATUS 1 MMG_VERSION_STATUS_1)

  IF ( MMG_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CI_DIR}/mmg.version MMG_MD5)

    IF ( NOT (${OLD_MMG_MD5} MATCHES ${MMG_MD5}) )
      SET ( GET_MMG_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMG_VERSION_STATUS_1}.
 Try to get it at the following link:
        https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmRktsVkFDTGlfdzQ
 then untar it in the ${CI_DIR} directory.")
  ENDIF()

  #--------------> mmg2d
  IF ( EXISTS ${CI_DIR}/mmg2d.version )
    FILE(MD5 ${CI_DIR}/mmg2d.version OLD_MMG2D_MD5)
  ELSE ( )
    SET ( OLD_MMG2D_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=1Lnvh7AldwEXS7WRa1VxsRqI7Xu7CgJNj
    ${CI_DIR}/mmg2d.version
    STATUS MMG2D_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMG2D_VERSION_STATUS 0 MMG2D_VERSION_STATUS_0)
  LIST(GET MMG2D_VERSION_STATUS 1 MMG2D_VERSION_STATUS_1)

  IF ( MMG2D_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CI_DIR}/mmg2d.version MMG2D_MD5)

    IF ( NOT (${OLD_MMG2D_MD5} MATCHES ${MMG2D_MD5}) )
      SET ( GET_MMG2D_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMG2D_VERSION_STATUS_1}.
 Try to get it at the following link:
        https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmV3BlUER4M0Z4MGs
 then untar it in the ${CI_DIR} directory.")
  ENDIF()

    #--------------> mmgs
  IF ( EXISTS ${CI_DIR}/mmgs.version )
    FILE(MD5 ${CI_DIR}/mmgs.version OLD_MMGS_MD5)
  ELSE ( )
    SET ( OLD_MMGS_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmSWtGV295a28wU2c
    ${CI_DIR}/mmgs.version
    STATUS MMGS_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMGS_VERSION_STATUS 0 MMGS_VERSION_STATUS_0)
  LIST(GET MMGS_VERSION_STATUS 1 MMGS_VERSION_STATUS_1)

  IF ( MMGS_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CI_DIR}/mmgs.version MMGS_MD5)

    IF ( NOT (${OLD_MMGS_MD5} MATCHES ${MMGS_MD5}) )
      SET ( GET_MMGS_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMGS_VERSION_STATUS_1}.
 Try to get it at the following link:
       https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmSWtGV295a28wU2c
 then untar it in the ${CI_DIR} directory.")
  ENDIF()

  #--------------> mmg3d
  IF ( EXISTS ${CI_DIR}/mmg3d.version )
    FILE(MD5 ${CI_DIR}/mmg3d.version OLD_MMG3D_MD5)
  ELSE ( )
    SET ( OLD_MMG3D_MD5 "0" )
  ENDIF ( )

  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmSWtGV295a28wU2c
    ${CI_DIR}/mmg3d.version
    STATUS MMG3D_VERSION_STATUS
    INACTIVITY_TIMEOUT 5)
  LIST(GET MMG3D_VERSION_STATUS 0 MMG3D_VERSION_STATUS_0)
  LIST(GET MMG3D_VERSION_STATUS 1 MMG3D_VERSION_STATUS_1)

  IF ( MMG3D_VERSION_STATUS_0 MATCHES 0)
    FILE(MD5 ${CI_DIR}/mmg3d.version MMG3D_MD5)

    IF (NOT (${OLD_MMG3D_MD5} MATCHES ${MMG3D_MD5}))
      SET ( GET_MMG3D_TESTS "TRUE" )
    ENDIF ()
  ELSE( )
    MESSAGE(WARNING "Failed to load a simple text file, download status:"
      " ${MMG3D_VERSION_STATUS_1}.
 Try to get it at the following link:
        https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmSWtGV295a28wU2c
 then untar it in the ${CI_DIR} directory.")
  ENDIF()

ENDIF()

# Download the tests if needed
#--------------> mmg
IF ( GET_MMG_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmg test cases. May be very long...")
  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=1Kd2aow6nfBI1i5dSN6lXMxaDKLrtpd6r
    ${CI_DIR}/mmg.tgz
    SHOW_PROGRESS)
  IF ( NOT EXISTS ${CI_DIR}/mmg.tgz )
    MESSAGE("\n")
    MESSAGE(WARNING "Fail to automatically download the mmg test cases.
Try to get it at the following link:
       https://drive.google.com/uc?export=download&id=1Kd2aow6nfBI1i5dSN6lXMxaDKLrtpd6r then untar it in the ${CI_DIR} directory.")
  ELSE()
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar xzf
      ${CI_DIR}/mmg.tgz
      WORKING_DIRECTORY ${CI_DIR}/
      )
    IF ( NOT EXISTS ${CI_DIR}/mmg.tgz )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmg "
        "test cases directory (mmg.tgz).
Try to untar it by hand in the ${CI_DIR} directory.")
    ENDIF()
    FILE(REMOVE ${CI_DIR}/mmg.tgz)
  ENDIF ()
ENDIF ()


#--------------> mmg2d
IF ( GET_MMG2D_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmg2d test cases. May be very long...")
  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=1Lnvh7AldwEXS7WRa1VxsRqI7Xu7CgJNj
    ${CI_DIR}/mmg2d.tgz
    SHOW_PROGRESS)
  IF ( NOT EXISTS ${CI_DIR}/mmg2d.tgz )
    MESSAGE("\n")
    MESSAGE(WARNING "Fail to automatically download the mmg2d test cases.
Try to get it at the following link:
       https://drive.google.com/uc?export=download&id=1Lnvh7AldwEXS7WRa1VxsRqI7Xu7CgJNj
then untar it in the ${CI_DIR} directory.")
  ELSE()
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar xzf
      ${CI_DIR}/mmg2d.tgz
      WORKING_DIRECTORY ${CI_DIR}/
      )
    IF ( NOT EXISTS ${CI_DIR}/mmg2d.tgz )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmg2d "
        "test cases directory (mmg2d.tgz).
Try to untar it by hand in the ${CI_DIR} directory.")
    ENDIF()
    FILE(REMOVE ${CI_DIR}/mmg2d.tgz)
  ENDIF ()
ENDIF ()

#--------------> mmgs
IF ( GET_MMGS_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmgs test cases. May be very long...")
  FILE(DOWNLOAD https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmcVdZb1EzaTR3ZlU
    ${CI_DIR}/mmgs.tgz
    SHOW_PROGRESS)
  IF ( NOT EXISTS ${CI_DIR}/mmgs.tgz )
    MESSAGE("\n")
    MESSAGE(WARNING "Fail to automatically download the mmgs test cases.
Try to get it at the following link:
       https://drive.google.com/uc?export=download&id=0B3X6EwOEKqHmcVdZb1EzaTR3ZlU
then untar it in the ${CI_DIR} directory.")
  ELSE()
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar xzf
      ${CI_DIR}/mmgs.tgz
      WORKING_DIRECTORY ${CI_DIR}/
      )
    IF ( NOT EXISTS ${CI_DIR}/mmgs.tgz )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmgs "
        "test cases directory (mmgs.tgz).
Try to untar it by hand in the ${CI_DIR} directory.")
    ENDIF()
    FILE(REMOVE ${CI_DIR}/mmgs.tgz)
  ENDIF ()
ENDIF ()

#--------------> mmg3d
SET(ADDRESS
  https://drive.google.com/uc?export=download&id=1WJK8mbFh81QFsDuOUt7kcavJwLO4yatO
  https://drive.google.com/uc?export=download&id=1SvznS9n57f1jIVoeFMM7-WVGJUu3OdMz
  https://drive.google.com/uc?export=download&id=1wACP1jut6Dz4mTf6uQARW7Koc1Zs_O9H
  https://drive.google.com/uc?export=download&id=142BueykwzDGS_Ne_RzJfgMGqQNhzJb6a
  https://drive.google.com/uc?export=download&id=1dEXKIApQiEkplI03bgVbThKMKTL4P0M_
  https://drive.google.com/uc?export=download&id=1KA5H7oS9HrtXT3sUpGU78YU7BwRHMVk7
  https://drive.google.com/uc?export=download&id=1duHPrEjdHrb1k9VwoV_-uYw4CPbMN9AM
  https://drive.google.com/uc?export=download&id=179k-asjM88ewVumZSQ9eUMWLJRfbQwIz
  https://drive.google.com/uc?export=download&id=1yjGvVah-vFNhwsImrHA0Bu5sIo1Fo5wW
  https://drive.google.com/uc?export=download&id=1PpQpC0OvJUTieTs0jd_A0Qb5EVP3QfZw
  https://drive.google.com/uc?export=download&id=1DbI0CCIYvDX-cPLZehazMV5wmu-K9K4b
  https://drive.google.com/uc?export=download&id=1MDALbXmXSpoVaHo4ghU6QpmQhsBP8Pjy
  https://drive.google.com/uc?export=download&id=1iil5UBwVgpcErUcKd_wJi_Oguuuyy0HT
  )

SET(FILENAME
  ${CI_DIR}/mmg3d.tgz.aa
  ${CI_DIR}/mmg3d.tgz.ab
  ${CI_DIR}/mmg3d.tgz.ac
  ${CI_DIR}/mmg3d.tgz.ad
  ${CI_DIR}/mmg3d.tgz.ae
  ${CI_DIR}/mmg3d.tgz.af
  ${CI_DIR}/mmg3d.tgz.ag
  ${CI_DIR}/mmg3d.tgz.ah
  ${CI_DIR}/mmg3d.tgz.ai
  ${CI_DIR}/mmg3d.tgz.aj
  ${CI_DIR}/mmg3d.tgz.ak
  ${CI_DIR}/mmg3d.tgz.al
  ${CI_DIR}/mmg3d.tgz.am
  )


IF ( GET_MMG3D_TESTS MATCHES "TRUE" )
  MESSAGE("-- Download the mmg3d test cases. May be very long...")

  SET(LOAD_OK 1)

  FOREACH( i RANGE 0 12)
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
then untar it in the ${CI_DIR} directory.")
      SET ( LOAD_OK 0 )
      BREAK()
    ENDIF()

  ENDFOREACH()

  IF ( ${LOAD_OK} )
    EXECUTE_PROCESS(
      COMMAND cat ${FILENAME}
      COMMAND tar -xzf -
      WORKING_DIRECTORY ${CI_DIR}/
      TIMEOUT 10000
      )
    IF ( NOT EXISTS ${CI_DIR}/mmg3d )
      MESSAGE("\n")
      MESSAGE(WARNING "Fail to automatically untar the mmg3d"
        "test cases directory (mmg3d.tgz.*).
Try to untar it by hand in the ${CI_DIR} directory: "
        "cat mmg3d.tgz.* | tar xzvf - ")
    ENDIF()

    FILE(REMOVE ${FILENAME})
  ENDIF()

ENDIF()
