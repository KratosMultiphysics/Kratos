# Install script for directory: /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX")
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/pastix" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/include/pastix/config.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/libpastix.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/include/pastix.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/pastix" TYPE FILE FILES
    "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/include/pastix/api.h"
    "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/include/pastix/old_api.h"
    "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/include/pastix/datatypes.h"
    "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/include/pastix/nompi.h"
    "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/include/pastix/order.h"
    "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/graph/graph.h"
    "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/include/cblas.h"
    "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/include/lapacke.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/lib/pkgconfig/pastix.pc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/lib/pkgconfig/pastixf.pc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/bin/pastix_env.sh")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/cmake_install.cmake")
  include("/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/cmake_install.cmake")
  include("/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/refinement/cmake_install.cmake")
  include("/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/cmake_install.cmake")
  include("/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/docs/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
