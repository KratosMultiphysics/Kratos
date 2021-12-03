## Shallow water application

This is a research application that provides a set of tools for oceanographic and hydrographic simulations over shallow domains.

### Implemented models

|               | BDF                | Crank-Nicolson     | Adams-Moulton      |
|---------------|:------------------:|:------------------:|:------------------:|
| Gravity waves | :heavy_check_mark: | :heavy_check_mark: |                    |
| Boussinesq    | :heavy_check_mark: |                    | :heavy_check_mark: |
| Shallow water | :heavy_check_mark: |                    |                    |

### Dependencies

This application does not have dependencies at compile time.

If the coupling with a Navier-Stokes is required, add the following applications to compilation:
- PfemFluidDynamicsApplication
- HDF5Application
- MappingApplication

### References

- M. Masó, I. De-Pouplana, E. Oñate. A FIC-FEM stabilized formulation for the shallow water equations over partially dry domains. Computer Methods in Applied Mechanics and Engineering, 389C (2022) 114362.
