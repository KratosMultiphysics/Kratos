## Metis Application

The Metis application provides an interface to the [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) graph partitioning library. This is used within Kratos Multiphysics to generate mesh partitions for MPI runs.

The applicaiton is organized in a series of processes, each one of them providing an interface to a different partitioning algorithm. The most commonly used one for finite elements is ```metis_divide_heterogeneous_input_process.h```, which generates a partition of the mesh based on the nodal graph and can be used for meshes combining multiple element types.

### Compilation
From Ubuntu 18.04 onwards, Metis can be installed with the following command:

```Shell
sudo apt-get install libmetis-dev
```

If not automatically detected, use

`-DMETIS_ROOT_DIR=String`

to specify the root directory for Metis library.
