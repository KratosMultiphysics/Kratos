# Install script for directory: /home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/lib/libmmg2d.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmg2d" TYPE FILE FILES
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg2d/libmmg2d.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg2d/libmmg2df.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypes.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypesf.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/bin/mmg2d_O3")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/lib/libmmgs.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmgs" TYPE FILE FILES
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmgs/libmmgs.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmgs/libmmgsf.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypes.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypesf.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/bin/mmgs_O3")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/lib/libmmg3d.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmg3d" TYPE FILE FILES
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg3d/libmmg3d.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg3d/libmmg3df.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypes.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypesf.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/bin/mmg3d_O3")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/lib/libmmg.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmg2d" TYPE FILE FILES
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg2d/libmmg2d.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg2d/libmmg2df.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmgs" TYPE FILE FILES
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmgs/libmmgs.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmgs/libmmgsf.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmg3d" TYPE FILE FILES
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg3d/libmmg3d.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg3d/libmmg3df.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg" TYPE FILE FILES
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypes.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypesf.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg/libmmg.h"
    "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg/libmmgf.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/vicente/bin/Kratos/applications/MeshingApplication/custom_external_libraries/mmg/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
