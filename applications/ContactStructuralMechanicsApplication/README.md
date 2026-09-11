# Contact Structural Mechanics Application

 |             **Application**             |                                                                                    **Description**                                                                                    |                              **Status**                              | **Authors** |
|:---------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------:|:-----------:|
| `ContactStructuralMechanicsApplication` | The *Contact Structural Mechanics Application* contains the contact mechanics implementations that can be used by the *Structural Mechanics Application* and *Constitutive Laws Application* within *Kratos Multiphysics* | <img src="https://img.shields.io/badge/Status-%F0%9F%94%A7Maintained-blue"  width="300px"> | [*Vicente Mataix Ferrándiz*](mailto:vicente.mataix-ferrandiz@siemens.com)  <br /> [*Alejandro Cornejo Velázquez*](mailto:acornejo@cimne.upc.edu)  |

<p align="center">
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/validation/double_arch/data/result.gif" alt="Double arch, frictionless" style="width: 300px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/validation/double_arch/data/result_frictional.gif" alt="Double arch, frictional" style="width: 300px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/use_cases/in_ring/data/animation.gif" alt="Cylinder in ring" style="width: 300px;"/>
 <img src="https://github.com/KratosMultiphysics/Examples/raw/master/contact_structural_mechanics/use_cases/hyperelastic_tubes/data/half_cylinders.gif" alt="Hyperelastic tubes" style="width: 300px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/use_cases/cylinders/data/horizontal_movement_2_frictional.gif" alt="Contacting cylinders, frictional" style="width: 300px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/use_cases/self_contact/data/animation.gif" alt="Self contact" style="width: 300px;"/>
</p>

The application implements **mortar (segment-to-segment) contact and mesh tying with dual Lagrange multipliers** for the structural solvers of Kratos: frictionless and frictional contact enforced with an augmented Lagrangian, a penalty or a multipoint-constraint approach, solved with a semi-smooth Newton method with a consistent (automatically generated) linearisation. The formulation is the one developed in the PhD thesis of Vicente Mataix Ferrándiz (UPC, 2020), Chapter 4.

## 😎 Features

- **Mesh tying conditions based on the mortar formulation** (dual Lagrange multipliers, generic tied variable, optional static condensation)
- **Augmented Lagrangian contact conditions based on the mortar formulation**
    * *Frictionless formulation* with a scalar multiplier (contact pressure) or a vector multiplier (*components* formulation, condensable)
    * *Frictional formulation* (Coulomb law, objective slip increment, pure-slip option)
    * *Axisymmetric* 2D variants
    * *Normal variation* (`NV`) variants with the derivatives of the normals in the tangent
- **Penalty contact conditions based on the mortar formulation**
    * *Frictionless* and *frictional* formulations, implicit and **explicit** dynamics
- **Simplified MPC conditions based on the mortar formulation**: the mortar weights build multipoint constraints, giving a simplified *NTN* / *NTS* behaviour
    * *Frictionless* and *frictional* formulations
    * *Mesh tying formulation with tension check*
- **Consistent linearisation** of the mortar operators, dual shape functions and normals, with the tangent matrices **generated symbolically** (sympy)
- **Self-contact compatible** (automatic master/slave assignment)
- **Contact search** with KD-tree / octree broad phase, oriented bounding boxes (separating axis theorem), mortar-consistent gap, dynamic search
- **Strategies, processes, builder-and-solvers, convergence criteria and a mixed U–LM linear solver** (static condensation of the dual multipliers) used by the contact formulation
- **Adaptive remeshing** strategies for contact (Hessian and SPR-error metrics with MMG, requires the `MeshingApplication`)
- **Frictional laws** (**WIP**: Coulomb wired, Tresca available as a class)
- **+120 Python tests** (small / nightly / validation, including validation benchmarks) **and 91 C++ tests**

## 🧭 Formulations at a glance

| Formulation | `contact_settings.mortar_type` | Python process | Extra DoFs | Friction | 2D / 3D / axisym |
|---|---|---|---|---|---|
| ALM, scalar multiplier | `ALMContactFrictionless` | `alm_contact_process` (`Frictionless`) | 1 per slave node | – | 2D, 3D, axisym |
| ALM, vector multiplier | `ALMContactFrictionlessComponents` | `alm_contact_process` (`FrictionlessComponents`) | dim per slave node | – | 2D, 3D |
| ALM frictional | `ALMContactFrictional[PureSlip]` | `alm_contact_process` (`Frictional[PureSlip]`) | dim per slave node | Coulomb | 2D, 3D, axisym |
| Penalty | `PenaltyContactFrictionless`, `PenaltyContactFrictional[PureSlip]` | `penalty_contact_process`, `explicit_penalty_contact_process` | none | Coulomb | 2D, 3D, axisym |
| MPC contact | – (`mpc_contact_settings`) | `mpc_contact_process` | none (constraints) | Coulomb | 2D, 3D |
| Mesh tying | `ScalarMeshTying`, `ComponentsMeshTying` | `mesh_tying_process` | tied components per slave node | – | 2D, 3D |

![Contact formulations at a glance](https://raw.githubusercontent.com/KratosMultiphysics/Kratos/master/docs/pages/Applications/Contact_Structural_Mechanics_Application/General/images/csma_formulation_matrix.png)

## 🚀 Quick start

Compile Kratos with the application (it depends on the `StructuralMechanicsApplication`; the `ConstitutiveLawsApplication` is optional for non-linear materials and the `MeshingApplication` for adaptive remeshing):

```sh
export KRATOS_APPLICATIONS="${KRATOS_APP_DIR}/StructuralMechanicsApplication;${KRATOS_APP_DIR}/ContactStructuralMechanicsApplication"
```

A contact simulation is a standard structural `ProjectParameters.json` plus two blocks: `contact_settings` inside `solver_settings` (which makes the structural solver wrapper pick the contact solver) and a contact process in `processes.contact_process_list`:

```json
"solver_settings" : {
    "solver_type"           : "Static",
    "analysis_type"         : "non_linear",
    "contact_settings"      : { "mortar_type" : "ALMContactFrictionless" },
    "convergence_criterion" : "contact_residual_criterion"
},
"processes" : {
    "contact_process_list" : [{
        "python_module" : "alm_contact_process",
        "kratos_module" : "KratosMultiphysics.ContactStructuralMechanicsApplication",
        "process_name"  : "ALMContactProcess",
        "Parameters"    : {
            "model_part_name"     : "Structure",
            "contact_model_part"  : { "0" : ["Contact_Part_1", "Contact_Part_2"] },
            "assume_master_slave" : { "0" : ["Parts_Parts_Auto2"] },
            "contact_type"        : "Frictionless"
        }
    }]
}
```

```python
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

with open("ProjectParameters.json") as f:
    parameters = KratosMultiphysics.Parameters(f.read())
StructuralMechanicsAnalysis(KratosMultiphysics.Model(), parameters).Run()
```

`contact_model_part` lists the sub-model-parts holding the skin conditions of the potentially contacting surfaces; `assume_master_slave` names the master side (leave it empty for self-contact). Results worth printing: `AUGMENTED_NORMAL_CONTACT_PRESSURE` (effective contact pressure), `WEIGHTED_GAP`, `LAGRANGE_MULTIPLIER_CONTACT_PRESSURE` (scaled by `SCALE_FACTOR`), the flags `ACTIVE`, `SLAVE`, `MASTER`. A complete worked example (2D Hertz benchmark) is in the [tutorial](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Tutorial_Hertz_2D.html).

## 📚 Documentation

The full documentation lives in the Kratos documentation site (sources in [`docs/pages/Applications/Contact_Structural_Mechanics_Application/`](../../docs/pages/Applications/Contact_Structural_Mechanics_Application)):

| Section | Pages |
|---|---|
| General | [Overview](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/General/Overview.html) · [Getting started](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/General/Getting_Started.html) |
| Theory | [Contact problem and state of the art](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Contact_Problem_And_State_Of_The_Art.html) · [Constrained optimisation methods](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Constrained_Optimisation_Methods.html) · [Frictionless contact](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Frictionless_Contact.html) · [Frictional contact](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Frictional_Contact.html) · [Mortar integration and dual Lagrange multipliers](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Mortar_Integration_And_Dual_Lagrange_Multipliers.html) · [Mesh tying](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Mesh_Tying.html) · [Linearisation and derivatives](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Linearisation_And_Derivatives.html) · [Automatic differentiation](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Automatic_Differentiation.html) |
| Contact search | [Search pipeline and bounding volumes](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Contact_Search/Search_Pipeline_And_Bounding_Volumes.html) · [Gap computation](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Contact_Search/Gap_Computation.html) · [Self contact](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Contact_Search/Self_Contact.html) |
| Implementation | [Architecture](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Architecture.html) · [Conditions](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Conditions.html) · [Strategies and convergence criteria](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Strategies_And_Convergence_Criteria.html) · [Builder-and-solvers and linear solvers](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Builder_And_Solvers_And_Linear_Solvers.html) · [Processes](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Processes.html) · [Utilities](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Utilities.html) · [Frictional laws and MPC constraint](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Frictional_Laws_And_MPC_Constraint.html) · [Variables and flags](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Variables_And_Flags_Reference.html) |
| Usage | [Solver settings](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Solver_Settings_Reference.html) · [Contact process settings](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Contact_Process_Settings_Reference.html) · [Tutorial (Hertz 2D)](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Tutorial_Hertz_2D.html) · [Output and post-processing](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Output_And_Postprocessing.html) · [Tips, troubleshooting and limitations](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Tips_Troubleshooting_And_Limitations.html) |
| Validation | [Benchmarks](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Validation/Benchmarks.html) · [Test suite](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Validation/Test_Suite_Reference.html) |
| Examples | [Applications gallery](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Examples/Applications_Gallery.html) · [Adaptive remeshing](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Examples/Adaptive_Remeshing.html) |
| Reference | [Bibliography](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Reference/Bibliography.html) · [Glossary](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Reference/Glossary.html) |

Each source folder also has a short `README.md` describing its content: [`custom_conditions`](custom_conditions/README.md), [`custom_strategies`](custom_strategies/README.md), [`custom_processes`](custom_processes/README.md), [`custom_utilities`](custom_utilities/README.md), [`custom_frictional_laws`](custom_frictional_laws/README.md), [`custom_linear_solvers`](custom_linear_solvers/README.md), [`custom_master_slave_constraints`](custom_master_slave_constraints/README.md), [`custom_python`](custom_python/README.md), [`python_scripts`](python_scripts/README.md), [`automatic_differentiation`](automatic_differentiation/README.md), [`tests`](tests/README.md).

The theory is developed in Chapter 4 of the *PhD thesis* authored by [Vicente Mataix Ferrándiz](mailto:vmataix@altair.com), available on [UPC Commons](https://upcommons.upc.edu/bitstream/2117/328952/1/TVMF1de1.pdf); the documentation reproduces its figures and equations and maps them to the code. The Doxygen documentation of the C++ classes is generated from `documents/doxyfile`.

## 🗂 Structure of the application

![Architecture of the application](https://raw.githubusercontent.com/KratosMultiphysics/Kratos/master/docs/pages/Applications/Contact_Structural_Mechanics_Application/General/images/csma_architecture_layers.png)

| Folder | Content |
|---|---|
| `custom_conditions/` | `PairedCondition`, `MortarContactCondition` and the ALM / penalty / axisymmetric families (AD-generated), `MeshTyingMortarCondition`, `MPCMortarContactCondition` — 68 registered conditions |
| `custom_strategies/` | Newton–Raphson / line-search / MPC contact strategies, 18 convergence criteria (active set, displacement/LM residuals, mesh error), 3 builder-and-solvers (header-only) |
| `custom_linear_solvers/` | `MixedULMLinearSolver`: static condensation of the dual multipliers |
| `custom_processes/` | Contact search (KD-tree / octree + OBB, simple / advanced activation, MPC variant), normal gap and normal check, ALM initialisation, penalty adaptation, dynamic factor, SPR error |
| `custom_utilities/` | Contact, active-set, self-contact, derivative, explicit-contribution and interface utilities |
| `custom_frictional_laws/` | `FrictionalLaw`, Coulomb, Tresca (WIP) |
| `custom_master_slave_constraints/` | `ContactMasterSlaveConstraint` (MPC route) |
| `custom_python/` | pybind11 bindings and `ProcessFactoryUtility` |
| `python_scripts/` | Contact solvers (static, implicit, explicit, MPC, adaptive remeshing), contact processes (ALM, penalty, explicit penalty, MPC, mesh tying), criteria factory, sympy helpers |
| `automatic_differentiation/` | sympy generators, C++ templates and theory notes of the generated conditions (sympy 1.2) |
| `tests/` | Python test suites (small / nightly / validation) and C++ gtests |
| `documents/` | Doxygen configuration |

## ⚙️ Examples

Examples can be found [here](https://github.com/KratosMultiphysics/Examples/tree/master/contact_structural_mechanics) (validation benchmarks: double arch, Hertz, full Hertz, shallow ironing, press fit; use cases: contacting cylinders, ironing with die, cylinder in ring, hyperelastic tubes, tooth model, arc pressing block, gears with plasticity, self-contact). They are described, together with the thesis benchmarks, in the [Applications gallery](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Examples/Applications_Gallery.html) and [Benchmarks](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Validation/Benchmarks.html) pages.

## 🧪 Tests

```sh
cd applications/ContactStructuralMechanicsApplication/tests
python3 test_ContactStructuralMechanicsApplication.py -l small      # also: nightly, validation, all
python3 <kratos>/kratos/python_scripts/testing/run_cpp_tests.py      # C++ tests (KRATOS_BUILD_TESTING=ON)
```

See [`tests/README.md`](tests/README.md) and the [test suite reference](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Validation/Test_Suite_Reference.html).

## 📈 Development timeline

![Development timeline](https://raw.githubusercontent.com/KratosMultiphysics/Kratos/master/docs/pages/Applications/Contact_Structural_Mechanics_Application/General/images/csma_timeline.png)

Created in August 2016 as a mortar contact prototype; the frictional, penalty, self-contact and remeshing capabilities were added in 2019, the MPC route in 2020; since 2022 the application is in maintenance mode.

## 📝 How to cite

If you use this application, please cite the thesis and Kratos:

```bibtex
@phdthesis{MataixFerrandiz2020,
  author = {Mataix Ferr{\'a}ndiz, Vicente},
  title  = {Innovative mathematical and numerical models for studying the deformation of shells during industrial forming processes with the Finite Element Method},
  school = {Universitat Polit{\`e}cnica de Catalunya},
  year   = {2020},
  url    = {https://upcommons.upc.edu/handle/2117/328952}
}

@article{Dadvand2010,
  author  = {Dadvand, Pooyan and Rossi, Riccardo and O{\~n}ate, Eugenio},
  title   = {An object-oriented environment for developing finite element codes for multi-disciplinary applications},
  journal = {Archives of Computational Methods in Engineering},
  volume  = {17},
  pages   = {253--297},
  year    = {2010},
  doi     = {10.1007/s11831-010-9045-2}
}
```

The mortar formulation with dual Lagrange multipliers follows A. Popp, *Mortar Methods for Computational Contact Mechanics and General Interface Problems*, PhD thesis, TU München, 2012.
