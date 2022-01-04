## Shallow water application

This is a research application that provides a set of tools for oceanographic and hydrographic simulations over shallow domains.

### Implemented models

|               | BDF                | Crank-Nicolson     | Adams-Moulton      |
|---------------|:------------------:|:------------------:|:------------------:|
| Gravity waves | :heavy_check_mark: | :heavy_check_mark: |                    |
| Boussinesq    | :heavy_check_mark: |                    | :heavy_check_mark: |
| Saint-Venant  | :heavy_check_mark: |                    |                    |

### Dependencies

This application does not have other application dependencies at compile time.

If the coupling with a Navier-Stokes is required, add the following applications to compilation:
- PfemFluidDynamicsApplication
- HDF5Application
- MappingApplication

The following Python libraries may be required:
- `numpy` is used by the wave generator, the benchmarks and other auxiliary processes
- `scipy` is used by the benchmarks
- `mpmath` is used by the solitary wave benchmark
- `h5py` is used in the convergence output
- `pandas` is used in some post process tools


### References

- M. Masó, I. De-Pouplana, E. Oñate. A FIC-FEM stabilized formulation for the shallow water equations over partially dry domains. Computer Methods in Applied Mechanics and Engineering, 389C (2022) 114362.
