## Poromechanics Application

The Poromechanics Application contains developments in coupled solid-pore fluid interaction problems within Kratos Multiphysics.

### Features:

- UPw small displacement element for saturated porous media (with
equal order interpolation, unstable under incompressible-undrained
conditions)

- Stable UPw small displacement element for saturated porous media
(with higher order interpolation for displacements)

- FIC-Stabilized UPw small displacement element for saturated porous media
(with equal order interpolation for displacements)

- UPw Quasi-zero-thickness interface elements for defining cracks and
joints

- Local linear elastic damage model (Simo-Ju and modified Von Mises)

- Non-local linear elastic damage model (Simo-Ju and modified Von
Mises)

- Bilinear cohesive fracture model (for quasi-zero-thickness interface elements)

- Fracture propagation utility based on the combination of the
damage model with the insertion of interface elements after remeshing
with GiD


### How to use MPI in Poromechanics Application

Make sure that the following lines are properly set in the configure.sh file:

> -DMETIS_APPLICATION=ON \\
>
> -DMETIS_INCLUDE_DIR="/usr/include/" \\
>
> -DUSE_METIS_5=ON \\
>
> -DPARMETIS_ROOT_DIR="/usr/lib/" \\
>
> -DTRILINOS_APPLICATION=ON \\
>
> -DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu/" \\
>
> -DTRILINOS_INCLUDE_DIR="/usr/include/trilinos/" \\
>
> -DTRILINOS_LIBRARY_PREFIX="trilinos_" \\
>
> -DEXTERNAL_SOLVERS_APPLICATION=ON \\
>
> -DCONSTITUTIVE_MODELS_APPLICATION=ON \\
>
> -DSOLID_MECHANICS_APPLICATION=ON \\
>
> -DFLUID_DYNAMICS_APPLICATION=ON \\
>
> -DPOROMECHANICS_APPLICATION=ON \\
>
> -DUSE_PORO_MPI=ON \\

Uncomment the following line in
[KratosPoromechanics.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/PoromechanicsApplication/custom_problemtype/Poromechanics_Application.gid/KratosPoromechanics.py).

> import KratosMultiphysics.TrilinosApplication as TrilinosApplication

*Note*: For the moment, MPI only works in Linux and requires compiling METIS_APPLICATION and TRILINOS_APPLICATION. Non-local Damage and Fracture Propagation features do not work in MPI.
