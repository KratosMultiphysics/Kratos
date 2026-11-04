# MedApplication

The Med Application reads and writes med-files, which contain mesh, field results and other data. Med-files are [HDF5](https://www.hdfgroup.org/solutions/hdf5/) files following a documented group/dataset layout; this application talks to that layout directly through the HDF5 C API, without depending on the external MED-library. This format is used by [Salome](https://www.salome-platform.org/) and [Code_Aster](https://code-aster.org).

## Installation
Only HDF5 is required to compile the application, e.g. via `sudo apt-get install libhdf5-dev` on Ubuntu.

## Usage
- In Salome, mesh groups are translated into SubModelParts. Different geometries and nodes can be added.
- SubSub ... Modelparts can be created by specifying a name with `.`. I.e. like it usually works in Kratos
- The number of characters is restricted in Med: 64 for main mesh name, and 80 for groups. Everything beyond these limits is cut.


## Development

- Use [HDFView](https://www.hdfgroup.org/downloads/hdfview/) or `h5dump` to inspect the med-files.
- Make sure to check the return value of every HDF5 function call.
- HDF5 does not check whether the data written to the file conforms to the MED layout. This must be ensured by the code writing it; see `custom_io/med_model_part_io.cpp` and `custom_io/med_hdf5_utilities.h` for the group/dataset layout MED expects.
