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

### Residual convergence with master-slave constraints

For the C++ GeoMechanics solving-strategy factory, the residual criterion selects
the residual representation from `strategy_type`:

- `quasi_newton` uses a constraint-reduced residual. Master DOFs remain active in
  the convergence norm and slave DOFs are excluded.
- Other strategies use a raw, untransformed residual by default. Both master and
  slave DOFs are excluded from the convergence norm.

The `residual_representation` solver setting can override this default:

```json
{
    "strategy_type": "newton_raphson",
    "convergence_criterion": "residual_criterion",
    "residual_representation": "constraint_reduced"
}
```

Supported values are `raw` and `constraint_reduced`. Choose the value that matches
the residual vector actually supplied to the criterion; this setting only selects
which DOFs contribute to its norm and does not transform the residual or reduce
the system. It applies to residual-based criteria, including those nested in
`and_criterion` and `or_criterion`.

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


