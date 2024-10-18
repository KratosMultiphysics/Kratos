# Shallow water application

This is a research application that provides a set of tools for oceanographic and hydrographic simulations over shallow domains. The background of the stabilization method is explained in [^1].

## Overview

|               | BDF                | Crank-Nicolson     | Adams-Moulton      |
|---------------|:------------------:|:------------------:|:------------------:|
| Gravity waves | :heavy_check_mark: | :heavy_check_mark: |                    |
| Saint-Venant  | :heavy_check_mark: |                    |                    |
| Boussinesq    | :heavy_check_mark: |                    | :heavy_check_mark: |

## Dependencies

This application does not have other application dependencies at compile time. The following Python libraries may be required:

- `scipy` is used by the wave generator and by the benchmarks
- `numpy` is used to generate solitary waves and analytical solutions by the benchmarks

If the coupling with the Navier-Stokes equations is required [^2], add the following applications to compilation:

- [HDF5Application](../HDF5Application/README.md)
- [MappingApplication](../MappingApplication/README.md)
- [PfemFluidDynamicsApplication](../PfemFluidDynamicsApplication/README.md)

## References

[^1]: M. Masó, I. De-Pouplana, E. Oñate. A FIC-FEM stabilized formulation for the shallow water equations over partially dry domains. Computer Methods in Applied Mechanics and Engineering, 389C (2022) 114362 [10.1016/j.cma.2021.114362](https://doi.org/10.1016/j.cma.2021.114362)

[^2]: M. Masó, A. Franci, I. de-Pouplana, A. Cornejo and E. Oñate, A Lagrangian-Eulerian procedure for the coupled solution of the Navier-Stokes and shallow water equations for landslide-generated waves. Advanced Modelling and Simulation in Engineering Sciences, (2022) [10.21203/rs.3.rs-1457837/v1](https://doi.org/10.21203/rs.3.rs-1457837/v1) (in press)
