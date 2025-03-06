# Install script for directory: /home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd

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
    set(CMAKE_INSTALL_CONFIG_NAME "")
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

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/opensubdiv/osd" TYPE FILE PERMISSIONS OWNER_READ GROUP_READ WORLD_READ FILES
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/bufferDescriptor.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/cpuEvaluator.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/cpuPatchTable.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/cpuVertexBuffer.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/mesh.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/nonCopyable.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/opengl.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/types.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/patchBasisTypes.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/patchBasis.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/cpuGLVertexBuffer.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/glLegacyGregoryPatchTable.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/glPatchTable.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/glVertexBuffer.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/glMesh.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/glslPatchShaderSource.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/glXFBEvaluator.h"
    "/home/basti/KratosMultiphysics/OpenSubdiv/opensubdiv/osd/glComputeEvaluator.h"
    )
endif()

