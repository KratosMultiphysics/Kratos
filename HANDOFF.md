# HANDOFF — Structural response extraction and upcoming Python refactoring

Prepared 2026-11-05 (Europe/Berlin). Repository root: `/home/alanany/repos/Kratos`.
Working application: `applications/StructuralMechanicsApplication`.

**Current instruction: create this handoff only. Do not refactor or modify Python source yet.**
The next refactoring's detailed scope has not been specified. Obtain that scope before
starting implementation; do not treat the preceding binding proposal as approval to implement it.

Evidence labels used below:
- **Verified current state:** inspected in the working tree while preparing this document.
- **Previously verified:** inspected or exercised during the conversation; rerun relevant checks
  before relying on it after changes. Historical in-memory tests were not saved as test files.
- **Decision:** explicit user requirement or a convention communicated to and accepted by the user.
- **Proposal/assumption:** not implemented or not independently established; keep distinct from facts.

## 1. Project purpose, structure, and conventions

Kratos Multiphysics is a C++20 finite-element/multiphysics framework with Python interfaces
through pybind11, built with CMake. `kratos/` contains the core, including entities, geometry,
quadrature, and Python bindings. `applications/` provides discipline-specific extensions.
StructuralMechanicsApplication contains structural elements, constitutive laws, conditions,
solvers, Python helpers, and tests.

Relevant layout, relative to the repository root:

```text
AGENTS.md
agents/instructions/{python-conventions,cpp-conventions,testing,cmake}.instructions.md
kratos/{includes,geometries,integration,python,utilities,tests}/
applications/StructuralMechanicsApplication/
    custom_elements/solid_elements/
    custom_constitutive/
    custom_python/
    python_scripts/helper/
    python_scripts/structural_components/panel/panel_response_extractor.py
    tests/test_responses/
    tests/test_helper_responses.py
    tests/test_StructuralMechanicsApplication.py
```

Read root `AGENTS.md` first and the matching `agents/instructions/` files before editing.
Keep changes surgical, retain neighboring style, use snake_case for new Python functions,
and use `KratosMultiphysics.KratosUnittest.TestCase` for Python tests. Never modify vendored
`external_libraries/` without explicit authorization. Test the most specific behavior first.

## 2. Target file and current responsibilities

**Verified current state:** the target is
`applications/StructuralMechanicsApplication/python_scripts/helper/responses.py`.
It contains `ResponseExtractor` and eight public methods, all with `mode="average"`:

| Public method | Per-element result |
| --- | --- |
| `extract_von_mises` | Nonnegative von Mises stress |
| `extract_signed_von_mises` | Von Mises magnitude with dominant-principal-stress sign |
| `extract_tresca` | Largest principal-stress difference |
| `extract_maximum_shear_stress` | Half the Tresca stress |
| `extract_signed_maximum_shear_stress` | Signed maximum shear stress |
| `extract_principle_stresses` | `[sigma1, sigma2, sigma3]`, descending algebraic order |
| `extract_principle_strains` | `[epsilon1, epsilon2, epsilon3]`, descending algebraic order |
| `extract_von_mises_strain` | Generic Abaqus/CAE MISES tensor invariant applied to strain |

**Decision:** preserve the existing public spelling `principle` in these two method names;
renaming to `principal` during refactoring would change the API. Private helper names use
`principal`. Scalar extractors return NumPy arrays; principal extractors return lists of
three-value lists. The outer order is `model_part.Elements` iteration order, not a dictionary
keyed by element ID. Empty results are an empty NumPy array or an empty list respectively.
The helper post-processes the current element state and does not run an analysis itself.

The shared `_extract_response` loop validates modes, identifies the continuum analysis,
selects a center index before requesting response data when necessary, requests GP output,
handles restricted plane-strain stress recovery, and reduces GP results. It accepts scalar
or vector-valued response callbacks, a signed-magnitude maximum flag, and a strain flag.
Other helpers classify analyses, normalize tensor/Voigt formats, compute eigenvalues and
invariants, and read the center-GP catalogue.

A support-summary comment block at the very beginning documents APIs, supported elements,
formats, modes, and restrictions. **Decision: update that block whenever support changes.**

## 3. Important decisions and extraction modes

**Decisions/current implementation:** mode strings are case-sensitive and exactly:

- `average`: compute the selected response at each GP, then take its arithmetic mean.
  Principal values are sorted at each GP and averaged by rank. Signed scalar values may
  cancel. Do not replace this with an invariant/eigendecomposition of an averaged tensor.
- `Gp_max`: unsigned scalar responses use the largest scalar; signed responses select the
  largest absolute response and retain that GP's sign. Equal magnitudes select the first GP.
  Principal stresses select the complete triplet at the GP with the largest **sigma1**;
  principal strains analogously use the largest **epsilon1**. Do not independently maximize
  all three components or mix values from different GPs. The principal-stress choice was
  explicitly confirmed by the user and must be clear in the opening comments.
- `center_gp`: use the stored zero-based index for the actual standard integration rule.
  A known rule without a center warns through `KM.Logger.PrintWarning`, then calls the
  `center_element` placeholder. Unknown element/rule combinations raise an unsupported-case
  error. A center response's returned GP-data length is checked against the active GP count.
- `center_element`: intentionally unimplemented; raises `NotImplementedError` when reached
  for a nonempty supported model part. Do not silently interpolate, average, or select a
  nearest point. Invalid mode strings raise `ValueError` before iteration.

### Center-GP catalogue and preserved old implementation

**Verified current state:** `python_scripts/helper/center_gauss_point_lookup.py` is a data-only
Python file derived from the matching `kratos_gauss_point_reference.md` and `.json`.
It contains 16 geometry types, 50 registered ordinary solid element names, six geometry
families, and 30 standard rules (`GI_GAUSS_1` through `GI_GAUSS_5` for each family).
It retains interpolation orders, default integration rules, active-GP-count keys, integration
orders, center-presence flags, local center coordinates, and zero-based center indices.

`responses.py` loads the sibling file once using `runpy.run_path` and `Path(__file__)`.
This supports direct `importlib` loading of the source helper by the manual runner.
The data is Python because `cmake_modules/install_function.cmake` currently installs `.py`
files from Python directories, not arbitrary adjacent JSON files. Do not casually replace
this with an import that works only from an installed package or with an uninstalled data file.

The active lookup identifies the geometry through
`element.GetGeometry().GetGeometryType().name.removeprefix("Kratos_")`, verifies the
registered element name against that geometry's listed names, and selects the rule using
`len(element.GetIntegrationPoints())`. Counts are unique within each catalogue family.
It does not compare coordinates or use a tolerance.

| Family | Reference center `(xi, eta, zeta)` | Rules with center: rule number -> index |
| --- | --- | --- |
| Triangle | `(1/3, 1/3, 0)` | `1 -> 0`, `3 -> 3` |
| Quadrilateral | `(0, 0, 0)` | `1 -> 0`, `3 -> 4`, `5 -> 12` |
| Tetrahedron | `(1/4, 1/4, 1/4)` | `1 -> 0` |
| Hexahedron | `(0, 0, 0)` | `1 -> 0`, `3 -> 13`, `5 -> 62` |
| Prism | `(1/3, 1/3, 1/2)` | None |
| Pyramid | `(0, 0, -1/2)` | `1 -> 0` |

The target is the **reference-element area/volume centroid mapped through the geometry**.
It need not equal `Geometry.Center()` (a nodal average) or the physical centroid of a
curved/distorted element. Triangle rule 3 has four GPs and center index 3, not the middle
array entry. Index 0 is valid; check `is None`, never truthiness.

**Decision:** the former coordinate/tolerance implementation of `_get_center_gp_index`
must remain commented out as the clearly labeled OLD IMPLEMENTATION. It was verified
exactly preserved, not deleted. Do not revive it as a fallback.

### Latest integration-method proposal — NOT implemented

The last discussion was a read-only proposal. The user asked about exposing
`GetIntegrationMethod()` without changing core behavior or breaking other code.
**No binding implementation has been authorized or made.**

**Verified in local source/previous live probe:**
- `kratos/includes/element.h` already has virtual `Element::GetIntegrationMethod() const`.
- `BaseSolidElement` overrides it to return initialized `mThisIntegrationMethod`.
- The getter is not bound on Python elements. The enum is already exposed as
  `KM.GeometryData.IntegrationMethod`.
- `geometry.GetDefaultIntegrationMethod()` may differ from the active element rule.
- `Properties[SMA.INTEGRATION_ORDER]` is consumed at initialization and can later become
  stale relative to the active rule. A live quad probe initialized with order 3 retained
  9 GPs after its property changed to 5; only reinitialization produced 25 GPs.
- `Element.GetIntegrationPoints()` internally calls the C++ active getter but constructs
  global coordinates for every GP. Although the Python lookup reads only its length,
  it still incurs that C++ coordinate-list construction.

**Proposed application-only option (preferred in the previous answer):** add a thin
module function `SMA.GetElementIntegrationMethod(element)` in
`custom_python/add_custom_utilities_to_python.cpp`, using a const `Element&` lambda that
returns `rElement.GetIntegrationMethod()`, with `py::arg("element")`, a docstring, and an
explicit `#include "includes/element.h"`. The existing application module already calls
`AddCustomUtilitiesToPython`. No new entity registration or CMake hook should be needed.
The return is the existing enum by value and the call preserves virtual dispatch.

Subsequently, if approved, `responses.py` and the compact data could use element-name +
actual-rule-name keys directly. A new standalone binding test was proposed; the existing
application test-suite entry must not be changed without specific permission. This option
requires rebuilding StructuralMechanicsApplication. Consider preserving the current
count-based route for compatibility with older application binaries if that is required.

**Alternative, also not approved:** add `.def("GetIntegrationMethod", &Element::GetIntegrationMethod, ... )`
in `kratos/python/add_mesh_to_python.cpp` to allow `element.GetIntegrationMethod()` directly.
That changes the core Python binding and requires a core-binding rebuild. Do not implement
it under an assumption that the earlier proposal was approval.

## 4. Mathematics, formats, assumptions, and invariants to preserve

At every GP, stress is the symmetric Cauchy tensor. Let its principal values satisfy
`sigma1 >= sigma2 >= sigma3`, let `I` be the 3x3 identity, and let `A:B = sum_ij A_ij B_ij`.

```text
s = sigma - trace(sigma)/3 * I
von_mises = sqrt(3/2 * (s:s))
tresca = sigma1 - sigma3
maximum_shear = (sigma1 - sigma3)/2

j = argmax_i abs(sigma_i)
signed_von_mises = sign(sigma_j) * von_mises
signed_maximum_shear = sign(sigma_j) * maximum_shear
```

**Decision:** the sign uses the largest-absolute principal stress, not the trace/pressure.
Exact opposite-sign magnitude ties choose the positive principal because eigenvalues are
sorted descending before `argmax(abs(...))`. Hydrostatic states have zero equivalent/shear
magnitude even if the dominant principal stress is nonzero. Signed shear is a tension/
compression reporting convention, not the directional sign of traction on a chosen plane.

Principal values are eigenvalues of symmetric tensors, calculated with `np.linalg.eigvalsh`
and reversed into descending algebraic order. Principal strains satisfy
`det(epsilon - lambda*I) = 0`. Eigenvectors/direction vectors are not returned.

Stress Voigt ordering is Kratos ordering, not Abaqus/Standard six-component ordering:
`[xx, yy, zz, xy, yz, xz]`; stress shear components are physical, with no division by two.
In 2D, supported formats include `[xx, yy, xy]`, `[xx, yy, zz, xy]`, and 2x2 tensors.
Full 3x3 tensors and six-component vectors preserve all provided components.

For reduced plane stress, omitted stress components are zero by the plane-stress condition.
For reduced plane strain, sigma_zz is generally nonzero. Recovery is strictly:

```text
sigma_zz = nu * (sigma_xx + sigma_yy)
```

This is allowed only for names starting `SmallDisplacementElement` or
`SmallDisplacementBbarElement`, with exactly `SMA.LinearElasticPlaneStrain2DLaw`, and no
nonzero initial stress or strain reported at GPs. Other missing-sigma_zz cases raise
`NotImplementedError`. Do not generalize this elastic relation to plasticity, anisotropy,
thermal/eigenstrain, or finite deformation without a separately justified implementation.

The strain extractors request **`KM.GREEN_LAGRANGE_STRAIN_VECTOR`**. In the specifically
supported small-displacement formulations this is computed with infinitesimal kinematics
(`B*u`, or the Bbar variant). The variable name does not imply that finite-strain Green-
Lagrange or logarithmic strain extraction is supported by these functions.

Strain Voigt order is `[exx, eyy, ezz, gamma_xy, gamma_yz, gamma_xz]`.
Convert engineering shear to tensor shear: `epsilon_ij = gamma_ij/2` for `i != j`.
Three- and four-component 2D variants are analogous to stress but contain engineering
shear. Full tensor off-diagonals already hold physical shear and must not be halved again.

```text
epsilon_dev = epsilon - trace(epsilon)/3 * I
implemented_mises_strain = sqrt(3/2 * (epsilon_dev:epsilon_dev))
                        = sqrt(((epsilon1-epsilon2)^2 +
                                (epsilon2-epsilon3)^2 +
                                (epsilon3-epsilon1)^2)/2)
```

**User decision and evidence qualification:** the user requested the Abaqus definition.
The chosen implementation applies the generic Abaqus/CAE `MISES` tensor operation to an
engineering-strain field after shear conversion. The docs explicitly specify the generic
Mises invariant and engineering-tensor handling, but no separate universal equation labeled
"von Mises total strain" was located. This is a documented interpretation/inference and
has not been numerically cross-validated against Abaqus. Do not overstate that evidence.
It differs by a factor of 1.5 from `sqrt(2/3 * epsilon_dev:epsilon_dev)` equivalent strain,
and it is not accumulated equivalent plastic strain `PEEQ`.

Reduced plane-strain strain output uses `epsilon_zz = 0` by kinematics. Reduced plane-stress
strain output also uses `epsilon_zz = 0`, specifically to match Abaqus/CAE **stored field
output**, despite generally nonzero physical thickness strain. This was researched and
communicated to the user. Do not replace it with a Poisson thickness-strain reconstruction
without clarifying that the response meaning changes. Full supplied out-of-plane strain
components are retained.

Relevant primary/code sources from the conversation:
- Ansys Motion Theory, sections 7.1.3.5.1, .9, .10: principal strains, maximum shear, signed shear:
  https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/motion_theory/motion_theory_finite_elements.html
- Abaqus generic stress invariants:
  https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMODRefMap/simamod-c-conventions.htm
- Abaqus `FieldOutput`, `MISES`, and `isEngineeringTensor`:
  https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKERRefMap/simaker-c-fieldoutputpyc.htm
- Abaqus output identifiers, Principal Value Output and Tensor Output (`E33=0` field convention):
  https://docs.software.vt.edu/abaqusv2025/English/SIMACAEOUTRefMap/simaout-c-stdoutputvar.htm
- `kratos/utilities/math_utils.h`, `StrainVectorToTensor`: Kratos engineering shear conversion.
- `custom_elements/solid_elements/small_displacement.cpp`: infinitesimal strain from `B*u`.
- `custom_elements/solid_elements/base_solid_element.cpp`: stress/strain GP output, initialization.
- Conventional alternative equivalent strain (for comparison, not the implemented equation):
  https://damask-multiphysics.org/dev/documentation/reference/processing_tools/post-processing.html#damask.mechanics.equivalent_strain_Mises

## 5. Dependencies and relationships

The target depends on NumPy, `KratosMultiphysics as KM`, StructuralMechanicsApplication as
`SMA`, `pathlib.Path`, `runpy.run_path`, and the sibling compact lookup file.

- `python_scripts/structural_components/panel/panel_response_extractor.py` was the initial
  structural/style reference. It averages tensor components for its composite response;
  that averaging behavior is not the current scalar-response averaging contract.
- `python_scripts/helper/kratos_gauss_point_reference.md` and `.json` are user-provided
  authoritative local rule references; retain them. Previously their recorded source hashes
  were verified against the checkout. Their center definitions differ from nodal averaging.
- `python_scripts/helper/center_gauss_point_lookup.py` is the runtime subset. It lists normal
  SmallDisplacement, Bbar where registered, TotalLagrangian, and UpdatedLagrangian names
  sharing catalogue geometries. Being listed for quadrature does not establish strain support.
- `python_scripts/helper/__init__.py` currently contains only a helper-package docstring.
- `python_scripts/helper/responses_temp.py` is a legacy/alternative implementation with
  `get_von_mises(center_gp, center_average, center_interpolated, model_part)`. It includes
  interpolation and shell behavior absent from the active class. It existed before this
  work; do not merge, delete, or treat it as the current specification without approval.
- `tests/test_helper_responses.py` imports that old `get_von_mises` name from `helper.responses`,
  which the current target does not export. This is a verified API mismatch. Leave the user's
  existing test file alone until permission/scope explicitly includes fixing it.
- `tests/test_responses/Main.py` directly loads the working-tree target with `importlib`, then
  runs standard `StructuralMechanicsAnalysis`; it does not rely on an installed helper copy.
- `CMakeLists.txt` already installs Python scripts using `kratos_python_install_directory`.
  The data-only Python catalogue participates in that mechanism without a new CMake edit.

## 6. Support boundaries, existing issues, and testing

**Verified support guards:** stress analyses require local and working dimensions to agree.
2D needs a law declaring `PLANE_STRESS_LAW` or `PLANE_STRAIN_LAW`; 3D needs a law with working
dimension 3 and strain size 6. Full Cauchy stress output is not globally restricted to elastic
analyses, but arbitrary nonlinear formulations have not all been tested. Trusses, beams,
shells, membranes, and axisymmetric analyses are outside the supported response scope.

The strain guard is narrower: names starting `SmallDisplacementElement` or
`SmallDisplacementBbarElement3D`, with exactly one of `LinearElastic3DLaw`,
`LinearElasticPlaneStress2DLaw`, `LinearElasticPlaneStrain2DLaw`. Finite-strain formulations,
other laws/subclasses, and 2D Bbar strain extraction are rejected. 2D Bbar uses four-component
laws, incompatible with the currently allowed built-in three-component 2D elastic laws.

Known limitations/gaps to preserve or address only with authorization:
- `center_element` remains unfinished; no-center `center_gp` therefore raises after warning.
- `center_gp` accepts only names/rules in the compact standard-solid catalogue.
- GP count identifies a rule only within the declared standard `GI_GAUSS_1..5` scope.
  Arbitrary C++-configured Lobatto/custom/element-owned rules could share a count and cannot
  be certified by that count alone; unknown counts/names are rejected, but matching custom
  counts are an inherent scope caveat. Do not advertise general custom-rule support.
- Current getter/count paths still construct coordinates internally; the faster direct enum
  route is a proposal, not existing functionality.
- The persistent legacy test targets a different API; current manual runner only exercises
  von Mises stress. Its CLI choices omit `Gp_max` and it does not select other extractors.
- The runner asserts reaction balance, not numerical equality of stress against a reference.
- A few per-method docstrings still say "geometry center"; the catalogue/implementation
  defines the mapped reference-element centroid. Preserve the actual current definition.

### Manual six-case assets

All under `applications/StructuralMechanicsApplication/tests/test_responses/`, with shared
`Main.py` and `README.md`; each case has `ProjectParameters.json`, `Materials.json`, `plate.mdpa`:

| Case | Geometry mesh | Element |
| --- | --- | --- |
| `plane_stress_q4` | 10x4 structured | SmallDisplacementElement2D4N |
| `plane_stress_q8` | 10x4 structured | SmallDisplacementElement2D8N |
| `plane_strain_q4` | 10x4 structured | SmallDisplacementElement2D4N |
| `plane_strain_q8` | 10x4 structured | SmallDisplacementElement2D8N |
| `solid_h8` | 10x4x2 structured | SmallDisplacementElement3D8N |
| `solid_h20` | 10x4x2 structured | SmallDisplacementElement3D20N |

Units: mm, N, MPa, tonne. Plate length100, width20, thickness2 mm; steel E210000 MPa,
nu0.3, density7.85e-9 tonne/mm^3. Fix all active displacements at x=0; apply +x tension at
x=100, with other boundaries free. 2D line load200 N/mm and 3D surface load100 MPa each
produce4000 N. Matching quadratic load conditions include midside nodes. One linear-static
load step; `move_mesh_flag=false`; core `skyline_lu_factorization` solver.

From the application directory, when a user authorizes running analyses:

```bash
PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 python3 tests/test_responses/Main.py
PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 python3 tests/test_responses/Main.py --case solid_h20 --mode center_gp
```

The default is all six cases, `average`, von Mises stress only. First-order quad/hex defaults
lack a center GP; second-order defaults include one. Per-case `results/` outputs contain
stress CSV, nodal CSV, and summary JSON and are ignored by the local `.gitignore`.

**Previously verified by this session:** input parsing, mesh connectivity/volumes, loads'
boundary sizes, constraints, initialization for all six models; in-memory `KratosUnittest`
checks for invariants, signs, rotations, all shear directions, Voigt/tensor equivalence,
plane strain/plane stress conventions, all modes, actual first/second-order elements,
unsupported cases, initial-state rejection for stress recovery, and existing API regression.
The most recent feature pass had 10 test methods; the catalogue pass had 8 methods including
80 geometry/rule combinations and 45 live initialization/rule combinations across all six
families. Those transient scripts were not persisted. No full suite/compiled rebuild or
Abaqus numerical comparison was performed by this assistant.

**Verified current disk evidence:** all six case result directories now exist with saved
`summary_average.json` files (they were not produced during this handoff task). Current
recorded downstream von Mises means and reaction-balance errors are:

| Case | Downstream mean (MPa) | Force-balance error (N) |
| --- | ---: | ---: |
| plane_stress_q4 | 99.9999993868 | 1.14e-11 |
| plane_stress_q8 | 99.9999993665 | 1.18e-11 |
| plane_strain_q4 | 88.8819654375 | 2.55e-11 |
| plane_strain_q8 | 88.8819435634 | 1.01e-10 |
| solid_h8 | 100.0000557911 | 5.84e-10 |
| solid_h20 | 100.0000017572 | -3.18e-9 |

These are saved results, not independently rerun in this task; their producing code revision
is not established. They are useful candidate baselines. Far-field ideal stress is100 MPa
for plane stress/free 3D contraction and100*sqrt(1-nu+nu^2)=88.8819 MPa for plane strain.
Clamp/corner effects make a global maximum-equals-reference assertion inappropriate.

Before a behavior-preserving refactor, establish reproducible baselines with focused checks
covering all eight APIs, outputs and ordering, average vs tensor averaging, signed magnitude
maxima and ties, principal-triplet selection, engineering shear conversion, Abaqus field
zero-zz convention, stress sigma_zz recovery limits, curved/high-order center lookup,
triangle rule3 index3, zero-index validity, overrides, no-center fallback, unknown rules,
and response-count mismatch. New persistent tests or edits outside the approved source
scope require the user's authorization. Use `PYTHONDONTWRITEBYTECODE=1` for read-only probes.

## 7. Current Git status and uncommitted work

**Verified snapshot before creating HANDOFF.md:**
- Branch: `smapp/responses`.
- HEAD: `a99085dc8f4` — `small bug fixes`.
- No tracked unstaged differences; no staged differences.
- 28 untracked files, all within StructuralMechanicsApplication. All helper/source/test assets
  below are untracked user work: ordinary `git diff` will not show their contents.
- `tests/test_StructuralMechanicsApplication.py` has no current modifications and does not
  import the untracked legacy helper test.
- Ignored result folders and Python caches may exist; preserve them.
- This new root `HANDOFF.md` is the sole additional repository file created for this task.

Repository-relative status captured at preparation time:

```text
?? applications/StructuralMechanicsApplication/python_scripts/helper/__init__.py
?? applications/StructuralMechanicsApplication/python_scripts/helper/center_gauss_point_lookup.py
?? applications/StructuralMechanicsApplication/python_scripts/helper/kratos_gauss_point_reference.json
?? applications/StructuralMechanicsApplication/python_scripts/helper/kratos_gauss_point_reference.md
?? applications/StructuralMechanicsApplication/python_scripts/helper/responses.py
?? applications/StructuralMechanicsApplication/python_scripts/helper/responses_temp.py
?? applications/StructuralMechanicsApplication/tests/test_helper_responses.py
?? applications/StructuralMechanicsApplication/tests/test_responses/.gitignore
?? applications/StructuralMechanicsApplication/tests/test_responses/Main.py
?? applications/StructuralMechanicsApplication/tests/test_responses/README.md
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_strain_q4/Materials.json
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_strain_q4/ProjectParameters.json
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_strain_q4/plate.mdpa
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_strain_q8/Materials.json
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_strain_q8/ProjectParameters.json
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_strain_q8/plate.mdpa
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_stress_q4/Materials.json
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_stress_q4/ProjectParameters.json
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_stress_q4/plate.mdpa
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_stress_q8/Materials.json
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_stress_q8/ProjectParameters.json
?? applications/StructuralMechanicsApplication/tests/test_responses/plane_stress_q8/plate.mdpa
?? applications/StructuralMechanicsApplication/tests/test_responses/solid_h20/Materials.json
?? applications/StructuralMechanicsApplication/tests/test_responses/solid_h20/ProjectParameters.json
?? applications/StructuralMechanicsApplication/tests/test_responses/solid_h20/plate.mdpa
?? applications/StructuralMechanicsApplication/tests/test_responses/solid_h8/Materials.json
?? applications/StructuralMechanicsApplication/tests/test_responses/solid_h8/ProjectParameters.json
?? applications/StructuralMechanicsApplication/tests/test_responses/solid_h8/plate.mdpa
```

Target source SHA256 at preparation:
`3742ac5208d3785cf00ccaa990674a78229cf1a28283b4693e17d7d4ab59e0f0`.

Do not run `git clean`, discard/reset changes, overwrite untracked files, stage, or commit
without authorization. Do not infer ownership solely from untracked status.

## 8. Instructions for the next Codex session

1. Read root `AGENTS.md`, this handoff, the relevant conventions, and current target/data
   files. Recheck Git status because the user may have made further edits.
2. Ask for/receive the actual refactoring requirements before implementing the refactor.
   The preceding integration-method exposure was explicitly proposal-only.
3. Preserve public method names, defaults, return types, iteration order, mode spelling,
   sign conventions, all equations, strain field conventions, support limits, and errors
   unless the user specifically authorizes a behavior change. Maintain opening comments.
4. Standing user constraint: normally modify only `responses.py`; ask for confirmation
   before editing other existing files. The user explicitly authorized creating the compact
   lookup and the manual test folder in earlier tasks, but that is not blanket permission
   to alter arbitrary files. Never change `tests/test_StructuralMechanicsApplication.py`
   without a new explicit user instruction: the user specifically forbade that edit.
5. Keep the old center-GP method commented and labeled. Retain both original reference
   documents and all pre-existing legacy files/tests, even when incompatible.
6. Keep Kratos core behavior and compatibility intact. Any binding implementation needs
   separate authorization; the proposed app-level const getter wrapper is the narrowest
   route. Do not monkey-patch the global Python Element class or infer active integration
   state solely from mutable properties/defaults.
7. For future work, validate focused behavior first; do not present initialization-only
   checks, saved result files, source-based inference, or one vendor's sign convention as
   proof of complete numerical equivalence with Abaqus.
8. Environment as of this task: bash; Python3.12 and a working Release Kratos build were
   previously available. Sandbox writable roots are the application directory and `/tmp`.
   The Git repository root lies outside the default writable application scope; creating
   this root document requires the filesystem approval mechanism. Do not work around that
   restriction by silently placing a requested root file in the application directory.
9. This handoff task permits only documentation creation. No Python/C++ source, lookup,
   input mesh/material/parameter, test, or existing result file has been edited in this turn.
