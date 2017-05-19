# Install script for directory: /home/rzorrilla/Kratos/applications

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

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/MeshingApplication/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/ExternalSolversApplication/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/StructuralMechanicsApplication/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/convection_diffusion_application/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/FluidDynamicsApplication/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/ALEapplication/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/FSIapplication/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/metis_application/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/trilinos_application/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/SolidMechanicsApplication/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/MappingApplication/cmake_install.cmake")
  include("/home/rzorrilla/Kratos/cmake_build_debug/applications/CompressiblePotentialFlowApplication/cmake_install.cmake")

endif()

