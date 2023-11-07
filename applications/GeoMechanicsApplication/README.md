## Geo-Mechanics Application

The Geo-Mechanics Application contains features needed for common geotechnical/geomechanical applications within Kratos Multiphysics.

### Features:
- K<sub>0</sub> procedure, Quasi-static, dynamic

- Staged analysis

- Automatic time stepping

- 2D (plane strain and axisymmetric) and 3D UPw small displacement element for saturated and partially saturated porous media (with
equal order interpolation, unstable under incompressible-undrained
conditions)

- 2D (plane strain and axisymmetric) and 3D  Stable UPw small displacement element for saturated and partially saturated porous media
(with higher order interpolation for displacements)

- 2D (plane strain and axisymmetric) and 3D FIC-Stabilized UPw small displacement element for saturated and partially saturated porous media
(with equal order interpolation for displacements)

- UPw Quasi-zero-thickness interface elements for defining cracks and
joints under saturated and partially saturated conditions

- UPw Updated-Lagrangian element for saturated and partially saturated porous media (with
equal order interpolation, unstable under incompressible-undrained
conditions)

- Stable UPw Updated-Lagrangian element for saturated and partially saturated porous media
(with higher order interpolation for displacements)

- 2D and 3D truss and cable elements

- 2D curved beam elemens with 3 nodes

- 1D, 2D and 3D steady-state and transient groundwater flow elements

- Loading User Defined Soil Models (UDSM) dll/so, written in PLAXIS forrmat

- Loading User Materials (UMAT) dll/so, written in ABAQUS forrmat

### How to compile Geo-Mechanics Application

Make sure that the following lines are properly set in the configuration file:

#### Windows:
~~~
CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\GeoMechanicsApplication;
~~~

#### Linux:
~~~
add_app ${KRATOS_APP_DIR}/LinearSolversApplication;
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication;
add_app ${KRATOS_APP_DIR}/GeoMechanicsApplication;
~~~

#### Note: 
- MPI has not been tested and does not work.

- The UMAT/UDSM constitutive models are not included in this repository. Some practical constitutive models can be found at https://soilmodels.com for instance.



