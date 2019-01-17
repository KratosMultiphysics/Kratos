# Install script for directory: /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/mod_files/")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/libpastixf.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fsimple" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fsimple")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fsimple"
         RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/fsimple")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fsimple" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fsimple")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fsimple"
         OLD_RPATH "/home/josep_maria/Kratos/libs:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fsimple")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/examples/fsimple.f90")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/flaplacian" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/flaplacian")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/flaplacian"
         RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/flaplacian")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/flaplacian" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/flaplacian")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/flaplacian"
         OLD_RPATH "/home/josep_maria/Kratos/libs:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/flaplacian")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/examples/flaplacian.f90")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fstep-by-step" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fstep-by-step")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fstep-by-step"
         RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/fstep-by-step")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fstep-by-step" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fstep-by-step")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fstep-by-step"
         OLD_RPATH "/home/josep_maria/Kratos/libs:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fstep-by-step")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/examples/fstep-by-step.f90")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fmultilap" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fmultilap")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fmultilap"
         RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/fmultilap")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fmultilap" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fmultilap")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fmultilap"
         OLD_RPATH "/home/josep_maria/Kratos/libs:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/fmultilap")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/examples/fmultilap.f90")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/examples/test_seq.in")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/examples/test_mt.in")
endif()

