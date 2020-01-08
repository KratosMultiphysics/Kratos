# Install script for directory: /home/jig/Kratos

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/jig/Kratos/bin/Release")
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
  message(STATUS "Deleting: /home/jig/Kratos/bin/Release/KratosMultiphysics")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(REMOVE_RECURSE "/home/jig/Kratos/bin/Release/KratosMultiphysics")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  message(STATUS "Deleting: /home/jig/Kratos/bin/Release/libs")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(REMOVE_RECURSE "/home/jig/Kratos/bin/Release/libs")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/KratosMultiphysics" TYPE FILE MESSAGE_NEVER FILES "/home/jig/Kratos/kratos/python_interface/__init__.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/KratosMultiphysics" TYPE FILE MESSAGE_NEVER FILES "/home/jig/Kratos/kratos/python_interface/kratos_globals.py")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/jig/Kratos/build/Release/external_libraries/gidpost/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/kratos/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/ConvectionDiffusionApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/DamApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/DEMApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/DemStructuresCouplingApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/ExternalSolversApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/FluidDynamicsApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/FreeSurfaceApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/MeshingApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/ParticleMechanicsApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/PfemFluidDynamicsApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/PoromechanicsApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/SwimmingDEMApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/applications/StructuralMechanicsApplication/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/external_libraries/triangle/cmake_install.cmake")
  include("/home/jig/Kratos/build/Release/kratos/run_kratos/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/jig/Kratos/build/Release/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
