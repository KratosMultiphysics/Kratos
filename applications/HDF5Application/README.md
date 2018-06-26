# HDF5Application

The *HDF5Application* enables the serialization of a model part with or without MPI using the [HDF5 library](https://support.hdfgroup.org/HDF5/). The model part is stored in and HDF5 file, which can be used for:

* Viewing a model part with [HDFVIEW](https://support.hdfgroup.org/products/java/hdfview/)
* Scientific visualization with tools supporting [XDMF](http://www.xdmf.org/index.php/Main_Page). Tools, which are known to work, include [ParaView 5.4](https://www.paraview.org/) and [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/).
* Checkpointing (under development).
* Re-partitioning (not supported yet).

## Installing HDF5 (minimum version 1.8)
The HDF5 C libraries are used with the *HDF5Application*. If Kratos is configured with MPI then the parallel HDF5 library must be installed. Otherwise, the serial HDF5 library is used.

### Serial HDF5 Library
#### Ubuntu
```
    sudo apt-get install libhdf5-dev
```
### Parallel HDF5 Library
#### Ubuntu
```
    sudo apt-get install libhdf5-openmpi-dev
```

## Build Instructions

1. Install serial or parallel HDF5 library

2. Configure Kratos to build the *HDF5Application*

    ```
    -DHDF5_APPLICATION=ON \
    ```

3. Build Kratos

## Installing h5py
This package is needed to use some of the python-files, e.g. `create_xdmf_file.py`
```
    sudo apt-get install python-h5py / python3-h5py
```

## Note
The minimum version for the GCC compiler is **4.9**. This is because earlier version don't fully support *regular expressions*.
