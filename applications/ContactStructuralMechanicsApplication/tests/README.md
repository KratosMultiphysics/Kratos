# Tests

![Test suite map](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Validation/images/csma_test_suite_map.svg)

## Running

```sh
# Python (from this folder); levels: small (default) ⊂ nightly ⊂ validation ⊂ all
python3 test_ContactStructuralMechanicsApplication.py -l small
python3 test_ContactStructuralMechanicsApplication.py -l nightly
# C++ gtests (application built with -DKRATOS_BUILD_TESTING=ON)
python3 <kratos>/kratos/python_scripts/testing/run_cpp_tests.py
```

About 15 tests require the `ConstitutiveLawsApplication` (detected at run time, `has_CL_application`) and the adaptive-remeshing tests require the `MeshingApplication`; they are skipped otherwise.

## Layout

| Entry | Content |
|---|---|
| `test_ContactStructuralMechanicsApplication.py` | Assembles the suites (40 small, +46 nightly, +37 validation active tests; 5 validation tests commented out: Hertz sphere axisymmetric ×2, ironing ×2, multi-layer). |
| `test_symbolic_generation.py` | Guard of the sympy code generation: regenerates the `2D2N` ALM frictionless condition with `automatic_differentiation/mortar_condition_generator.py` and checks structure and numerical equivalence with the committed file (skipped when sympy is not installed). |
| `contact_structural_mechanics_test_factory.py` | `ContactStructuralMechanicsTestFactory`: each test class names a `file_name`; the factory reads `<file_name>_parameters.json`, runs a `StructuralMechanicsAnalysis` and checks the results with `from_json_check_result_process` against `<file_name>_results.json`. `frictionless_by_components = True` switches the case to the vector-multiplier formulation (`Components*` classes). |
| `SmallTests.py`, `NightlyTests.py`, `ValidationTests.py` | 32 / 55 / 25 factory classes. |
| `test_check_normals_process.py`, `test_double_curvature_integration.py`, `test_dynamic_search.py`, `test_process_factory.py` | Stand-alone unit tests (normal check, exact double-curvature integration, dynamic search, process factory). |
| `ALM_frictionless_contact_test_2D/`, `ALM_frictionless_contact_test_3D/`, `ALM_frictional_contact_test_2D/`, `ALM_frictional_contact_test_3D/`, `penalty_frictionless_contact_test_2D/`, `penalty_frictionless_contact_test_3D/`, `penalty_frictional_contact_test_2D/`, `mpc_contact_tests/`, `mesh_tying_test/` | Case data: `.mdpa`, `*_parameters.json`, `*_materials.json`, `*_results.json`. |
| `auxiliary_files_for_python_unittest/` | Meshes for the stand-alone unit tests. |
| `cpp_tests/` | 91 gtest cases in `KratosContactStructuralMechanicsFastSuite`: derivatives (49), weighted gap (11), MixedULM solver (8), active set (5), integration (4), search (3), self-contact (3), mesh-tying condition (2), interface preprocess (2), AALM / ALM variables / parent elements / contact utilities (1 each). |

## Categories

Patch tests (ALM, penalty, MPC, mesh tying; 2D/3D, matching and non-matching, sloped, mixed tri/quad), Taylor patch test, Hertz (2D/3D, frictionless/frictional), frictional stick/slip/threshold tests, beams and plates (MPC), large-displacement patch tests, mesh-moving cases, self-contact, multi-layer, explicit penalty, dynamic search, exact integration, mortar mapping and adaptive remeshing. The mapping between the thesis benchmarks and the test cases is in [Benchmarks](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Validation/Benchmarks.md).

## Full documentation

- [Test suite reference](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Validation/Test_Suite_Reference.html) · [source](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Validation/Test_Suite_Reference.md)
- [Benchmarks](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Validation/Benchmarks.html)
