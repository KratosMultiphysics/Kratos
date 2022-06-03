## Geo-Mechanics Application

The Geo-Mechanics Application contains features needed for common geotechnical/geomechanical applications within Kratos Multiphysics.

### Features:
- K<sub>0</sub> procedure, Quasi-static, dynamic

- Staged analysis

- Automatic time stepping

- UPw small displacement element for saturated porous media (with
equal order interpolation, unstable under incompressible-undrained
conditions)

- Stable UPw small displacement element for saturated porous media
(with higher order interpolation for displacements)

- FIC-Stabilized UPw small displacement element for saturated porous media
(with equal order interpolation for displacements)

- UPw Quasi-zero-thickness interface elements for defining cracks and
joints

- UPw Updated-Lagrangian element for saturated porous media (with
equal order interpolation, unstable under incompressible-undrained
conditions)

- Stable UPw Updated-Lagrangian element for saturated porous media
(with higher order interpolation for displacements)

- Reading and using dll/so of User Defined Soil Models (UDSM) based on PLAXIS forrmat

- Reading and using dll/so of UMAT based on ABAQUS forrmat

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

*Note*: For the moment, MPI does not work.
