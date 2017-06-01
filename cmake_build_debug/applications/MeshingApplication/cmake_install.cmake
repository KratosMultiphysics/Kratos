# Install script for directory: /home/rzorrilla/Kratos/applications/MeshingApplication

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
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosMeshingApplication.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosMeshingApplication.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosMeshingApplication.so"
         RPATH "/home/rzorrilla/Kratos/libs:/home/rzorrilla/Compiled_libraries/boost_1_59_0/stage/lib:/usr/lib/python3.5/config-3.5m-x86_64-linux-gnu/libpython3.5m.so:/usr/lib/libblas/libblas.so.3.6.0:/usr/lib/lapack/liblapack.so.3.6.0")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/libs" TYPE SHARED_LIBRARY FILES "/home/rzorrilla/Kratos/cmake_build_debug/applications/MeshingApplication/KratosMeshingApplication.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosMeshingApplication.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosMeshingApplication.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosMeshingApplication.so"
         OLD_RPATH "/home/rzorrilla/Compiled_libraries/boost_1_59_0/stage/lib:/usr/lib/python3.5/config-3.5m-x86_64-linux-gnu/libpython3.5m.so:/usr/lib/libblas/libblas.so.3.6.0:/usr/lib/lapack/liblapack.so.3.6.0:/home/rzorrilla/Kratos/cmake_build_debug/kratos:/home/rzorrilla/Kratos/cmake_build_debug/external_libraries/zlib:"
         NEW_RPATH "/home/rzorrilla/Kratos/libs:/home/rzorrilla/Compiled_libraries/boost_1_59_0/stage/lib:/usr/lib/python3.5/config-3.5m-x86_64-linux-gnu/libpython3.5m.so:/usr/lib/libblas/libblas.so.3.6.0:/usr/lib/lapack/liblapack.so.3.6.0")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosMeshingApplication.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/KratosMultiphysics" TYPE FILE FILES "/home/rzorrilla/Kratos/applications/MeshingApplication/MeshingApplication.py")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/MeshingApplication/custom_external_libraries/triangle/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/MeshingApplication/custom_external_libraries/tetgen/tetgen1.4.3/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/MeshingApplication/custom_external_libraries/tetMeshOpt/cmake_install.cmake")

endif()

