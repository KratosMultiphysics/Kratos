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
