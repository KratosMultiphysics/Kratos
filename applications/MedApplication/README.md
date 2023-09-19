# MedApplication

The Med Application an interface to the MED-library. This library writes med-files, which contain mesh, field results and other data, and is based on [HDF5](https://www.hdfgroup.org/solutions/hdf5/). This format is used by [Salome](https://www.salome-platform.org/) and [Code_Aster](https://code-aster.org).

## Installation
The MED-library is an external library, which must be installed before the application can be compiled

On Ubuntu, it can be installed with `sudo apt-get install libmedc-dev`. This installs all required dependencies, including HDF5

The source code is available on the Salome website for a manual installation. In this case also HDF5 needs to be installed separately.

Use `MED_ROOT` to specify the path to the MED installation in the CMake of Kratos.
