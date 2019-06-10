## Metis Application

The Metis application provides an interface to the [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) graph partitioning library. This is used within Kratos Multiphysics to generate mesh partitions for MPI runs.

The applicaiton is organized in a series of processes, each one of them providing an interface to a different partitioning algorithm. The most commonly used one for finite elements is ```metis_divide_heterogeneous_input_process.h```, which generates a partition of the mesh based on the nodal graph and can be used for meshes combining multiple element types.
