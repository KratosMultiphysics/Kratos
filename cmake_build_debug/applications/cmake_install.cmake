# Install script for directory: /home/laurin/kratos_rep_nov17/Kratos/applications

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
  include("/home/laurin/kratos_rep_nov17/Kratos/cmake_build_debug/applications/ExternalSolversApplication/cmake_install.cmake")
  include("/home/laurin/kratos_rep_nov17/Kratos/cmake_build_debug/applications/SolidMechanicsApplication/cmake_install.cmake")
  include("/home/laurin/kratos_rep_nov17/Kratos/cmake_build_debug/applications/PfemApplication/cmake_install.cmake")
  include("/home/laurin/kratos_rep_nov17/Kratos/cmake_build_debug/applications/PfemSolidMechanicsApplication/cmake_install.cmake")
  include("/home/laurin/kratos_rep_nov17/Kratos/cmake_build_debug/applications/PfemFluidDynamicsApplication/cmake_install.cmake")
  include("/home/laurin/kratos_rep_nov17/Kratos/cmake_build_debug/applications/ContactMechanicsApplication/cmake_install.cmake")
  include("/home/laurin/kratos_rep_nov17/Kratos/cmake_build_debug/applications/ConstitutiveModelsApplication/cmake_install.cmake")
  include("/home/laurin/kratos_rep_nov17/Kratos/cmake_build_debug/applications/StructuralMechanicsApplication/cmake_install.cmake")
  include("/home/laurin/kratos_rep_nov17/Kratos/cmake_build_debug/applications/constitutive_laws_application/cmake_install.cmake")

endif()

