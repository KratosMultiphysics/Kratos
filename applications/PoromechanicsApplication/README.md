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

### How to use MPI in Poromechanics Application (Linux only)

First, install the following packages from a terminal:

- trilinos-dev
- libmetis-dev
- libmetis5
- libopenmpi-dev
- libscotch-dev
- libtrilinos-amesos-dev
- libtrilinos-aztecoo-dev
- libtrilinos-epetra-dev
- libtrilinos-epetraext-dev
- libtrilinos-ifpack-dev
- libtrilinos-ml-dev
- libtrilinos-teuchos-dev
- openmpi-bin

Then, make sure that the following applications are added in the configure.sh file:

> add_app ${KRATOS_APP_DIR}/ExternalSolversApplication;
>
> add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication;
>
> add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication;
>
> add_app ${KRATOS_APP_DIR}/PoromechanicsApplication;
>
> add_app ${KRATOS_APP_DIR}/MetisApplication;
>
> add_app ${KRATOS_APP_DIR}/TrilinosApplication;

And also that the following options are set:

> -DUSE_MPI=ON \\
>
> -DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu" \\
>
> -DTRILINOS_INCLUDE_DIR="/usr/include/trilinos" \\
>
> -DTRILINOS_LIBRARY_PREFIX="trilinos_" \\

Uncomment the following line in
[KratosPoromechanics.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/PoromechanicsApplication/custom_problemtype/Poromechanics_Application.gid/KratosPoromechanics.py).

> import KratosMultiphysics.TrilinosApplication

*Note*: For the moment, MPI only works in Linux. Non-local Damage and Fracture Propagation features do not work in MPI.