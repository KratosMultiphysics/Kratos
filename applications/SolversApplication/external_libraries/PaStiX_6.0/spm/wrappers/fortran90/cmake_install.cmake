# Install script for directory: /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/libspmf.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_driver" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_driver")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_driver"
         RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/spm_driver")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_driver" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_driver")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_driver"
         OLD_RPATH "/home/josep_maria/Kratos/libs:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_driver")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/examples/spm_driver.f90")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_user" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_user")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_user"
         RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE EXECUTABLE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/spm_user")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_user" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_user")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_user"
         OLD_RPATH "/home/josep_maria/Kratos/libs:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib:/usr/lib/x86_64-linux-gnu:/home/josep_maria/Kratos/libs")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/examples/spm_user")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/examples" TYPE FILE FILES "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/examples/spm_user.f90")
endif()

