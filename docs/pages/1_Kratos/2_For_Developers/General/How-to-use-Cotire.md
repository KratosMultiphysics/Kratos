---
title: Cotire
keywords: 
tags: [How-to-use-Cotire.md]
sidebar: kratos_for_developers
summary: 
---

# What is Cotire?

[Cotire](https://github.com/sakra/cotire) (compile time reducer) is a CMake module that speeds up the build process of CMake based build systems by fully automating techniques as precompiled header usage and single compilation unit builds for C and C++.

Depending on factors like hardware, compiler, the number of files in the target and the complexity
of the C/C++ code, the build process of targets that use a cotire generated precompiled header
will be sped up from 10 to 40 percent.

A unity build may be up to 90 percent faster than the one file at a time build of the original
target. Single compilation unit builds however are very unlikely to work without source code
modifications, because they [break][http://altdevblog.com/2011/08/14/the-evils-of-unity-builds/] the use of some C and C++ language features.

Generally, modern C++ code which makes heavy use of header-only libraries will profit the most from
cotiring.

# How to compile with Cotire?

Add `-DUSE_COTIRE=ON` to your configure file. This makes the compiler use precompiled headers.
Then `substitute make install -j4` with `make all_unity install/fast -j4` to speed up the compilation.
It will only work if your application was already adapted to Cotire through CMakeLists.txt.

# How do I add compatibly with Cotire with my app?

To use cotire in your CMake project, add the following include directive to the beginning of the `CMakeList.txt` of your application:

    include(cotire)

To speed the build process of a CMake library or executable target, just apply the `cotire`
function to the target:

    add_executable(MyExecutable ${MyExecutableSources})
    target_link_libraries(MyExecutable ${MyExecutableLibraries})
    cotire(MyExecutable)

Cotire looks at the properties of the target provided by CMake (e.g., target type, source files,
compile flags, preprocessor defines, include directories, ...) and sets up custom commands that
will generate a unity source file, a prefix header and a precompiled header at build time
specially tailored to the target.

For the generation of the prefix header, cotire will automatically choose headers used by the
target that are outside of the project directory and thus are likely to change infrequently.
The precompiled prefix header is then applied to the target to speed up the compilation process.

To use an existing manually maintained prefix header instead of the automatically generated one,
set the `COTIRE_CXX_PREFIX_HEADER_INIT` property before invoking cotire:

    set_target_properties(MyExecutable PROPERTIES COTIRE_CXX_PREFIX_HEADER_INIT "stdafx.h")
    cotire(MyExecutable)

As a side effect, cotire generates a new target named `MyExecutable_unity`, which lets you perform
a unity build for the original target. The unity target inherits all build settings from the
original target, including linked library dependencies.

For Makefile based generators you can then invoke a unity build that produces the same output as
the original target, but does so much faster by entering:

    $ make MyExecutable_unity

See the advanced usage section of the [cotire manual][https://github.com/sakra/cotire/blob/master/MANUAL.md] for information on how to configure the cotire process (e.g., how to make the unity build use all available processor
cores).

For example in the case of the StructuralMechanicsAplication we have:

```cmake
set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosStructuralMechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/MeshingApplication )

## generate variables with the sources
set( KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE
  ## MAIN FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/structural_mechanics_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/structural_mechanics_application_variables.cpp

  ....(the others cpp files)
)

## generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## generate variables with the sources
set( KRATOS_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE
  ## CUSTOM PYTHON
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/structural_mechanics_python_application.cpp
  ....(the others custom python)
)

add_library(KratosStructuralMechanicsCore SHARED ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosStructuralMechanicsCore PUBLIC KratosCore)
set_target_properties(KratosStructuralMechanicsCore PROPERTIES COMPILE_DEFINITIONS "STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosStructuralMechanicsApplication MODULE THIN_LTO ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosStructuralMechanicsApplication PUBLIC KratosStructuralMechanicsCore)
set_target_properties(KratosStructuralMechanicsApplication PROPERTIES PREFIX "")

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN "*.csv")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

if(USE_COTIRE MATCHES ON)
    cotire(KratosStructuralMechanicsCore)
    cotire(KratosStructuralMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosStructuralMechanicsCore DESTINATION libs )
install(TARGETS KratosStructuralMechanicsApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/StructuralMechanicsApplication.py" DESTINATION KratosMultiphysics )
```

# CCache?

You can even reduce more the compilation time with [CCache](https://ccache.samba.org/) (which also reduces compilation time).

TODO: FINISH THIS



