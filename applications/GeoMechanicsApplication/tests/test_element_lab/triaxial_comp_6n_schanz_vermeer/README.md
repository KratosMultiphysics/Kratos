# Triaxial Schanz Vermeer elasticity test with 6 node axisymmetric elements

This staged lab-element test reuses the axisymmetric 2D6N triaxial setup and applies the Schanz Vermeer stress dependant incremental elastic law.

The specimen starts from a uniform in-situ stress state $\sigma_{xx}' = \sigma_{yy}' = \sigma_{zz}' = -100$ [Pa]. During both stages the lateral pressure is kept at -100 [Pa].

## Stage 1 - Unloading

The top normal stress is reduced from -100 Pa to -50 Pa.

## Stage 2 - Reloading

The top normal stress is increased from -50 Pa back to -100 Pa.

With Poisson's ratio $\nu$ = 0.0 [-] and

```math
E = E^{ref} ( \frac{-sigma_3'}{p_{ref}} )^m
```

using $E^{ref}$ = 1.0e7 [Pa], $p_{ref}$ = 50 [Pa] and $m$ = 1 [-], the analytical axial stage strain for the stress change between -100 and -50 Pa is:

```math
\epsilon_{yy} = \frac{p_{ref}}{E^{ref}} ln(2.0)
```

which is used to validate the top-node stage displacements.
