MORSE CMake modules
====================

This project provides a collection of CMake modules that can be shared
among projects using CMake as build system.

For now it is mainly constituted of "Find" modules that help detecting
installed libraries on the system. These modules are located in

Get morse_cmake
---------------------

To use latest development states of morse_cmake, please clone the
master branch:

    git clone git@gitlab.inria.fr:solverstack/morse_cmake.git

Documentation
---------------------

See the file [morse_cmakefind_doc.org](modules/find/morse_cmakefind_doc.org).

Installation
---------------------

We recommend to use this project as a `git submodule` of your project.

    # Example if morse_cmake is defined as a git submodule in ./cmake_modules/
    git submodule add https://gitlab.inria.fr/solverstack/morse_cmake.git cmake_modules/morse_cmake

To use MORSE modules you have to add the path to the modules in your
CMake project and include the MorseInit module:

    # Example if Morse CMake modules are located in ./cmake_modules/morse_cmake/modules
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules" )
    # Include the init module
    include(MorseInit)

Testing
---------------------

See the file [README](modules/find/tests/README.md).

Get involved!
---------------------

### Mailing list

TODO

### Contributions

https://gitlab.inria.fr/solverstack/morse_cmake/blob/master/CONTRIBUTING.md

### Authors

The following people contributed to the development of morse_modules:
  * Cedric Castagnede
  * Mathieu Faverge, PI
  * Florent Pruvost, PI

If we forgot your name, please let us know that we can fix that mistake.

### Licence

https://gitlab.inria.fr/solverstack/morse_cmake/blob/master/LICENCE.txt
