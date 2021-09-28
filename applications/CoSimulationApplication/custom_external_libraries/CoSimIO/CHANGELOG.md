# Changelog

All important and notable changes in the _CoSimIO_ will be documented in this file.

## 2.0.0
- Changes in Connecting:
    - Now instead of specifying `connection_name` directly, `my_name` and `connect_to` are used for the call to `Connect`.
    - The `CoSimIO::Info` that is returned from `Connect` now contains the `connection_name` that needs to be used in subsequent calls to _CoSimIO_.
    - This change was done to make the connection more robust and to make the selection of primary/secondary partner clear. This is only used internally e.g. for who opens a port and who connects to it.
- Introduced `CoSimIO::ModelPart` for exchanging of meshes
    - Simplified version of [`Kratos::ModelPart`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/includes/model_part.h)
    - Simplifies and unifies the usage of `Import-/ExportMesh`
    - See the tutorials on how to use it:
        - [C++](docs/model_part/model_part_cpp.md)
        - [C](docs/model_part/model_part_c.md)
        - [Python](docs/model_part/model_part_python.md)
- FileCommunication:
    - By default now done in folder. This way leftovers from previous simulations can be easily deleted (done automatically).
    - working directory can be specified
    - stability of initial connection was significantly improved.
- Python interface: Data is no longer copied when going from Python to C++ and vice versa.
    - `Import-/ExportData` now uses `CoSimIO::DoubleVector` (small wrapper around `std::wrapper`)
    - `Import-/ExportMesh` now uses `CoSimIO::ModelPart`
- Continuous Integration:
    - Adding Python 3.9 (now has Python v3.5 - v3.9)
    - Adding CentOS 7 build with GCC 4.8.5
    - Enforcing C89 standard

- Many improvements and cleanups under the hood
