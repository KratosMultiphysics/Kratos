# MedApplication

The Med Application an interface to the MED-library. This library writes med-files, which contain mesh, field results and other data, and is based on [HDF5](https://www.hdfgroup.org/solutions/hdf5/). This format is used by [Salome](https://www.salome-platform.org/) and [Code_Aster](https://code-aster.org).

## Installation
The MED-library is an external library, which must be installed before the application can be compiled

### Ubuntu

On Ubuntu, it can be installed with `sudo apt-get install libmedc-dev`. This installs all required dependencies, including HDF5

The source code is available on the Salome website for a manual installation. In this case also HDF5 needs to be installed separately.

Use `MED_ROOT` to specify the path to the MED installation in the CMake of Kratos.

### Arch / Manjaro

Packages related to *Salome* and *MED* for arch-based distros can be installed from the [AUR](https://en.wikipedia.org/wiki/Arch_Linux#Arch_User_Repository_(AUR)). The MedApplication requires [med-serial](https://aur.archlinux.org/packages/med-serial) (for non-MPI builds) or [med-openmpi](https://archlinux.org/packages/extra/x86_64/med-openmpi/) (for MPI builds with OpenMPI).
```
sudo pacman -S med-serial med-openmpi
```

## Usage
- In Salome, mesh groups are translated into SubModelParts. Different geometries and nodes can be added.
- SubSub ... Modelparts can be created by specifying a name with `.`. I.e. like it usually works in Kratos
- The number of characters is restricted in Med: 64 for main mesh name, and 80 for groups. Everything beyond these limits is cut.


## Development
- Use [HDFView](https://www.hdfgroup.org/downloads/hdfview/) to inspect the med-files.
- Make sure to check the return value of every med-library function call.
- The med library does not check if wrong data is written to the file. This must be ensured by the user, the med-library is a thin wrapper around HDF.
