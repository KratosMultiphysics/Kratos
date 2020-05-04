# Install script for directory: /home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg

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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg2d.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg2d.so.5.4.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg2d.so.5"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg2d.so.5.4.3"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg2d.so.5"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg2d.so.5.4.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg2d.so.5"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg2d.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg2d.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg2d.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg2d.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg2d.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg2d.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg2d.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmg2d" TYPE FILE FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg2d/libmmg2d.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/mmg2d/libmmg2df.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypes.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/common/libmmgtypesf.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xapplix" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/bin/mmg2d_O3")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg2d_O3")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmgs.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmgs.so.5.4.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmgs.so.5"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmgs.so.5.4.3"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmgs.so.5"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmgs.so.5.4.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmgs.so.5"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmgs.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmgs.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmgs.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmgs.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmgs.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmgs.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmgs.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmgs" TYPE FILE FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/mmgs/libmmgs.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/mmgs/libmmgsf.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypes.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/common/libmmgtypesf.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xapplix" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/bin/mmgs_O3")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmgs_O3")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg3d.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg3d.so.5.4.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg3d.so.5"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg3d.so.5.4.3"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg3d.so.5"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg3d.so.5.4.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg3d.so.5"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg3d.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg3d.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg3d.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg3d.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg3d.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg3d.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg3d.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmg3d" TYPE FILE FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg3d/libmmg3d.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/mmg3d/libmmg3df.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypes.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/common/libmmgtypesf.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xapplix" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/bin/mmg3d_O3")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mmg3d_O3")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg.so.5.4.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg.so.5"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg.so.5.4.3"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg.so.5"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg.so.5.4.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg.so.5"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/lib/libmmg.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmmg.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmg2d" TYPE FILE FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg2d/libmmg2d.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/mmg2d/libmmg2df.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypes.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/common/libmmgtypesf.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmgs" TYPE FILE FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/mmgs/libmmgs.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/mmgs/libmmgsf.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypes.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/common/libmmgtypesf.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg/mmg3d" TYPE FILE FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg3d/libmmg3d.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/mmg3d/libmmg3df.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/common/libmmgtypes.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/src/common/libmmgtypesf.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mmg" TYPE FILE FILES
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg/libmmg.h"
    "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/src/mmg/libmmgf.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/mmg/MmgTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/mmg/MmgTargets.cmake"
         "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/CMakeFiles/Export/lib/cmake/mmg/MmgTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/mmg/MmgTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/mmg/MmgTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/mmg" TYPE FILE FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/CMakeFiles/Export/lib/cmake/mmg/MmgTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/mmg" TYPE FILE FILES "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/CMakeFiles/Export/lib/cmake/mmg/MmgTargets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/mhashemi/programming/Git/Kratos-working/applications/MeshingApplication/custom_external_libraries/mmg/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
