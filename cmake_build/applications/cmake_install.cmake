# Install script for directory: /home/bsaridar/Software/Kratos_FullDebug/Kratos/applications

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
    set(CMAKE_INSTALL_CONFIG_NAME "FullDebug")
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

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/bsaridar/Software/Kratos_FullDebug/Kratos/cmake_build/applications/ExternalSolversApplication/cmake_install.cmake")
  include("/home/bsaridar/Software/Kratos_FullDebug/Kratos/cmake_build/applications/MeshingApplication/cmake_install.cmake")
  include("/home/bsaridar/Software/Kratos_FullDebug/Kratos/cmake_build/applications/StructuralMechanicsApplication/cmake_install.cmake")
  include("/home/bsaridar/Software/Kratos_FullDebug/Kratos/cmake_build/applications/FluidDynamicsApplication/cmake_install.cmake")
  include("/home/bsaridar/Software/Kratos_FullDebug/Kratos/cmake_build/applications/MeshMovingApplication/cmake_install.cmake")
  include("/home/bsaridar/Software/Kratos_FullDebug/Kratos/cmake_build/applications/EmpireApplication/cmake_install.cmake")
  include("/home/bsaridar/Software/Kratos_FullDebug/Kratos/cmake_build/applications/MetisApplication/cmake_install.cmake")
  include("/home/bsaridar/Software/Kratos_FullDebug/Kratos/cmake_build/applications/MappingApplication/cmake_install.cmake")
  include("/home/bsaridar/Software/Kratos_FullDebug/Kratos/cmake_build/applications/CompressiblePotentialFlowApplication/cmake_install.cmake")
  include("/home/bsaridar/Software/Kratos_FullDebug/Kratos/cmake_build/applications/CoSimulationApplication/cmake_install.cmake")

endif()

