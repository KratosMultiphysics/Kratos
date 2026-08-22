# Triaxial unloading-reloading test with 6 noded elements

This staged lab-element test reuses the axisymmetric 2D6N triaxial setup and applies the E_ur incremental elastic law.

The specimen starts from a uniform in-situ stress state of -100 Pa in xx, yy and zz. During both stages the lateral pressure is kept at -100 Pa.

## Stage 1 - Unloading

The top normal stress is reduced from -100 Pa to -50 Pa.

## Stage 2 - Reloading

The top normal stress is increased from -50 Pa back to -100 Pa.

With nu = 0 and

E = E_ur^ref ((-sigma_3') / p_ref)^m

using E_ur^ref = 1.0e7 Pa, p_ref = 50 Pa and m = 1, the analytical axial stage strain for the stress change between -100 and -50 Pa is:

epsilon_yy = (p_ref / E_ur^ref) ln(2)

which is used to validate the top-node stage displacements.
