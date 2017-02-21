# Install script for directory: /home/gcasas/kratos/external_libraries/zlib

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
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
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
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libz.so.1.2.8"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libz.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libz.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/home/gcasas/kratos/libs:/home/gcasas/compiled_libs/boost_1_57_0/stage/lib:/usr/lib/x86_64-linux-gnu/libpython3.4m.so:/usr/lib/atlas-base/libf77blas.so.3.0:/usr/lib/atlas-base/libatlas.so.3.0:/usr/lib/lapack/liblapack.so.3.0")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/libs" TYPE SHARED_LIBRARY FILES
    "/home/gcasas/kratos/cmake_build/external_libraries/zlib/libz.so.1.2.8"
    "/home/gcasas/kratos/cmake_build/external_libraries/zlib/libz.so.1"
    "/home/gcasas/kratos/cmake_build/external_libraries/zlib/libz.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libz.so.1.2.8"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libz.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libz.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
           NEW_RPATH "/home/gcasas/kratos/libs:/home/gcasas/compiled_libs/boost_1_57_0/stage/lib:/usr/lib/x86_64-linux-gnu/libpython3.4m.so:/usr/lib/atlas-base/libf77blas.so.3.0:/usr/lib/atlas-base/libatlas.so.3.0:/usr/lib/lapack/liblapack.so.3.0")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/libs" TYPE STATIC_LIBRARY FILES "/home/gcasas/kratos/cmake_build/external_libraries/zlib/libz.a")
endif()

