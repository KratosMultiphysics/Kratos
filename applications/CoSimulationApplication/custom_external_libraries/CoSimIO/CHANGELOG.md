# Changelog

All important and notable changes in the _CoSimIO_ are documented in this file.

## 2.0.0
- Changes in Connecting:
    - Now instead of specifying `connection_name` directly, `my_name` and `connect_to` are used for the call to `Connect`.
    - The `CoSimIO::Info` that is returned from `Connect` now contains the `connection_name` that needs to be used in subsequent calls to _CoSimIO_.
    - This change was done to make the connection more robust and to make the selection of primary/secondary partner clear. This is only used internally e.g. for who opens a port and who connects to it.
- Introduced `CoSimIO::ModelPart` for exchanging of meshes
    - Simplified version of [`Kratos::ModelPart`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/includes/model_part.h)
    - Simplifies and unifies the usage of `Import-/ExportMesh`
    - See the tutorials on how to use it:
        - [C++](https://kratosmultiphysics.github.io/CoSimIO/model_part/model_part_cpp.html)
        - [C](https://kratosmultiphysics.github.io/CoSimIO/model_part/model_part_c.html)
        - [Python](https://kratosmultiphysics.github.io/CoSimIO/model_part/model_part_python.html)
- FileCommunication:
    - By default now done in folder. This way leftovers from previous simulations can be easily deleted (done automatically).
    - working directory can be specified
    - stability of initial connection was significantly improved.
- Python interface: Data is no longer copied when going from Python to C++ and vice versa.
    - `Import-/ExportData` now uses `CoSimIO::DoubleVector` (small wrapper around [`std::vector`](https://en.cppreference.com/w/cpp/container/vector))
    - `Import-/ExportMesh` now uses `CoSimIO::ModelPart`
- Continuous Integration:
    - Adding Python 3.9 (now has Python v3.5 - v3.9)
    - Adding CentOS 7 build with GCC 4.8.5
    - Enforcing C89 standard

- Many improvements and cleanups under the hood

## 2.0.1
- Bugfix in remote-controlled CoSimulation (now settings can be passed to the registered functions)
- For this the `Info` object can now hold `Info` objects itself, hence making it possible to build hierarchical structures.

## 3.0.0
- Extensive documentation was added: https://kratosmultiphysics.github.io/CoSimIO/
- Now the C++ version of CoSimIO is not header-only any more. Check [here](https://kratosmultiphysics.github.io/CoSimIO/build_options.html) for the available build options.
- Support for MPI-parellism is now available (for all supported languages, i.e. C++, C and Python), along with extensive documentation and testing
- Pyramid topologies are available now. See [here](https://github.com/KratosMultiphysics/CoSimIO/pull/271) for details.
- Communication:
    - Communication via pipes was added (currently only supported in Linux). See [here](https://kratosmultiphysics.github.io/CoSimIO/communication.html#pipe-based-communication) for the details.
    - Communication now uses serialization for exchanging data.
    - FileCommuncation now longer adds an index to each file name
    - FileCommuncation: Instead of using `rename` for file synchronization, it is now possible to use auxiliary files. See [here](https://github.com/KratosMultiphysics/CoSimIO/pull/254) for the details.
- Improvements and extensions to `CoSimIO::ModelPart`:
    - The `ModelPart` interface was extended to support ghost nodes. See [here](https://kratosmultiphysics.github.io/CoSimIO/model_part/model_part_cpp.html#interface-for-distributed-modelparts-mpi) for details
    - ModelPart uses `intrusive_ptr` for storing nodes and elements to save memory
- C-Interface: Some functions like `CoSimIO_CreateInfo` or `CoSimIO_CreateModelPart` internally allocate memory and return a pointer to the allocated object. Discarding the return value leads to memory leaks. Compiler warnings were added to warn the user. Check [here](https://github.com/KratosMultiphysics/CoSimIO/pull/181) for more details.
- All CMake macros of _CoSimIO_ now start with `CO_SIM_IO_`. This affects especially `BUILD_C` and `BUILD_PYTHON` which are changed to `CO_SIM_IO_BUILD_C` and `CO_SIM_IO_BUILD_PYTHON`.
- Several minor improvements to the CI (continuous integration)
- Internal errors now give much better error messages including detailed stacktraces

## 4.0.0
- (Interprocess) Communication
    - Socket based communication was added (using TCP network sockets with IPv4). Documentation is available [here](https://kratosmultiphysics.github.io/CoSimIO/communication.html#socket-based-communication).
    Due to its versatility it is the **new default** (previously it was file-based communication)
    - Experimental support for unix domain sockets was added (see [here](https://kratosmultiphysics.github.io/CoSimIO/communication.html#unix-domain-socket-based-communication))
    - PipeCommunication now supports large data (data larger than pipe buffer). Furthermore the buffer size can be configured with `buffer_size`, see the [documentation](https://kratosmultiphysics.github.io/CoSimIO/communication.html#pipe-based-communication).
    - Experimental support for communication based on `MPI_Ports` was added (see [here](https://kratosmultiphysics.github.io/CoSimIO/communication.html#mpi-based-communication)). It is available when connection is done with MPI and can be enabled with `CO_SIM_IO_BUILD_MPI_COMMUNICATION`
- Mesh container (`ModelPart`)
    - Creating entities in the ModelPart is now significantly faster. Especially when creating many entities the speedup is several orders of magnitude.
    - New interfaces are added to the ModelPart with which multiple entities can be created.
- Improved and extended documentation
- Improved synchronization during initial handshake to avoid deadlocks
- Printing timing information for communication is now unified for all communication methods (can be enabled with `print_timing`)
- other minor interal improvements and fixes
