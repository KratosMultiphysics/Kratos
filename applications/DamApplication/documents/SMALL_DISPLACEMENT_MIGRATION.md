# DamApplication small-displacement migration: final architecture

This document records the completed migration of the DamApplication solid-mechanics
kernel to the StructuralMechanicsApplication `SmallDisplacement` element, and the
permanent architecture that results. It is the authoritative maintenance reference
for the decisions, policies and intentional compatibility differences introduced
by Phases 5A-6D.

---

## 1. Purpose of the refactor

Historically DamApplication maintained its own duplicated small-displacement solid
mechanical kernel (`SolidElement`, `SmallDisplacementElement`,
`SmallDisplacementThermoMechanicElement`). This duplicated kernel was removed.
StructuralMechanicsApplication `SmallDisplacement` now owns the generic solid
mechanics, while all Dam-specific constitutive, smoothing, nonlocal-orchestration,
construction, acoustic and coupled-workflow behavior remains in DamApplication.

The permanent architecture is:

```
SMA SmallDisplacement  ->  parameter-aware Dam constitutive laws
```

There is no Dam-owned small-displacement element class anymore.

---

## 2. Final ownership model

### StructuralMechanicsApplication owns (for Dam solid mechanics)

* small-displacement kinematics and B-matrix strain computation;
* stiffness assembly and internal forces;
* body-force assembly;
* element lifecycle interface;
* integration-point constitutive dispatch;
* consistent and lumped mass (M1);
* damping contribution;
* `Create` / `Clone` / `Check`;
* the generic constitutive-law integration-point output interface.

### DamApplication owns

**Constitutive behavior**

* thermal linear laws (`ThermalLinearElastic3DLaw`, 2D plane-strain/plane-stress);
* nodal thermal laws (`ThermalLinearElastic*Nodal`);
* local thermal damage (`ThermalLocalDamage3DLaw` + Simo-Ju 3D/2D);
* nonlocal thermal damage (`ThermalNonlocalDamage3DLaw` + Simo-Ju / Modified-Mises 3D/2D);
* specialized thermal/mechanical output semantics (see §10);
* automatic restart Properties rebinding (see §13).

**Post-processing / orchestration**

* `DamNodalCauchyStressExtrapolationProcess`;
* the smoothing schemes and `custom_utilities/dam_nodal_stress_smoothing_utilities.hpp`;
* `DamNonlocalDamageUtilities` (custom_utilities);
* nonlocal activation/orchestration;
* construction/selfweight-specific workflows.

**Acoustic / coupled workflows**

* `WaveEquationElement`;
* acoustic conditions (`FreeSurface*`, `InfiniteDomain*`, `AddedMass*`);
* `UPCondition`;
* `DamPScheme`; `DamUPScheme`;
* P/U-P/T-U-P solver orchestration.

**Compatibility**

* the 22 historical element-name aliases;
* the frozen historical binary restart fixture and its regression test.

---

## 3. Historical aliases

All 22 historical Dam small-displacement registration strings are **permanent
compatibility aliases**. Each resolves to the StructuralMechanics `SmallDisplacement`
prototype for the matching geometry. Mechanical and thermo-mechanical families are
no longer distinct element classes: they are two compatibility names for the SAME
SMA implementation, and the thermal vs mechanical behavior is governed solely by the
constitutive law assigned in `Properties`.

| Historical string | Geometry | SMA prototype | Runtime type | GP count | Family | Mass | `.mdpa` | binary restart |
|---|---|---|---|---|---|---|---|---|
| SmallDisplacementSolidElement2D3N | T3 | mSmallDisplacementElement2D3N | SmallDisplacement | 1 | mech | M1 | OK | OK |
| SmallDisplacementSolidElement2D4N | Q4 | mSmallDisplacementElement2D4N | SmallDisplacement | 4 | mech | M1(≡) | OK | OK |
| SmallDisplacementSolidElement2D6N | T6 | mSmallDisplacementElement2D6N | SmallDisplacement | 3 | mech | M1 | OK | OK |
| SmallDisplacementSolidElement2D8N | Q8 | mSmallDisplacementElement2D8N | SmallDisplacement | 9 | mech | M1(≡) | OK | OK |
| SmallDisplacementSolidElement2D9N | Q9 | mSmallDisplacementElement2D9N | SmallDisplacement | 9 | mech | M1(≡) | OK | OK |
| SmallDisplacementSolidElement3D4N | T4 | mSmallDisplacementElement3D4N | SmallDisplacement | 1 | mech | M1 | OK | OK |
| SmallDisplacementSolidElement3D6N | Prism6 | mSmallDisplacementElement3D6N | SmallDisplacement | 6 | mech | M1(≡) | OK | OK |
| SmallDisplacementSolidElement3D8N | H8 | mSmallDisplacementElement3D8N | SmallDisplacement | 8 | mech | M1(≡) | OK | OK |
| SmallDisplacementSolidElement3D10N | T10 | mSmallDisplacementElement3D10N | SmallDisplacement | 4 | mech | M1 | OK | OK |
| SmallDisplacementSolidElement3D15N | Prism15 | mSmallDisplacementElement3D15N | SmallDisplacement | 12 | mech | M1(≡) | OK | OK |
| SmallDisplacementSolidElement3D20N | H20 | mSmallDisplacementElement3D20N | SmallDisplacement | 27 | mech | M1(≡) | OK | OK |
| SmallDisplacementSolidElement3D27N | H27 | mSmallDisplacementElement3D27N | SmallDisplacement | 27 | mech | M1(≡) | OK | OK |
| SmallDisplacementThermoMechanicElement2D3N | T3 | mSmallDisplacementElement2D3N | SmallDisplacement | 1 | thermo | M1 | OK | OK |
| SmallDisplacementThermoMechanicElement2D4N | Q4 | mSmallDisplacementElement2D4N | SmallDisplacement | 4 | thermo | M1(≡) | OK | OK |
| SmallDisplacementThermoMechanicElement2D6N | T6 | mSmallDisplacementElement2D6N | SmallDisplacement | 3 | thermo | M1 | OK | OK |
| SmallDisplacementThermoMechanicElement2D8N | Q8 | mSmallDisplacementElement2D8N | SmallDisplacement | 9 | thermo | M1(≡) | OK | OK |
| SmallDisplacementThermoMechanicElement2D9N | Q9 | mSmallDisplacementElement2D9N | SmallDisplacement | 9 | thermo | M1(≡) | OK | OK |
| SmallDisplacementThermoMechanicElement3D4N | T4 | mSmallDisplacementElement3D4N | SmallDisplacement | 1 | thermo | M1 | OK | OK |
| SmallDisplacementThermoMechanicElement3D8N | H8 | mSmallDisplacementElement3D8N | SmallDisplacement | 8 | thermo | M1(≡) | OK | OK |
| SmallDisplacementThermoMechanicElement3D10N | T10 | mSmallDisplacementElement3D10N | SmallDisplacement | 4 | thermo | M1 | OK | OK |
| SmallDisplacementThermoMechanicElement3D20N | H20 | mSmallDisplacementElement3D20N | SmallDisplacement | 27 | thermo | M1(≡) | OK | OK |
| SmallDisplacementThermoMechanicElement3D27N | H27 | mSmallDisplacementElement3D27N | SmallDisplacement | 27 | thermo | M1(≡) | OK | OK |

`M1` = SMA default exact/elevated consistent mass (intentional for simplex);
`M1(≡)` = coincides with the historical Dam default for non-simplex geometries.

---

## 4. Mass policy M1

`M1 = StructuralMechanics default consistent-mass behavior`.

* **Quad / Hex / Prism**: the previous Dam default integration already produced the
  same consistent mass as the SMA exact/elevated rule, so behavior is numerically
  equivalent.
* **Simplex (T3, T6, T4, T10)**: the historical no-flag Dam behavior used
  under-integrated/default simplex mass; the aliases now use the SMA exact/elevated
  consistent mass. This is an **intentional numerical change** (M1).
* Lumped mass remains compatible.
* The obsolete Dam `COMPUTE_CONSISTENT_MASS_MATRIX` API was removed and no longer
  exists.

**Dynamic consequence**: the simplex consistent-mass change can affect inertia,
mass-proportional Rayleigh damping, modal/frequency behavior and transient
dynamics. Static/quasi-static analyses are unaffected when mass is unused;
non-simplex consistent mass was already equivalent; lumped mass is unchanged.

---

## 5. Restart policy R4

* **`.mdpa`**: historical Dam names are preserved permanently and resolve to SMA
  prototypes.
* **New binary restarts**: the runtime type is SMA; the Serializer writes the
  canonical SMA registration name (first-registration semantics).
* **Old production-format binary Dam restarts**: the frozen genuine legacy fixture
  (`tests/cpp_tests/fixtures/legacy_thermo_3d8n_damage.dat`) proves
  `historical archive -> historical stored name -> current production alias ->
  SMA runtime -> preserved constitutive state -> successful first post-load
  material response`.
* **ASCII trace mode**: cross-implementation compatibility is **not** guaranteed
  (the old hierarchy had an additional `BaseClass` nesting level). This is a
  documented caveat, not a supported compatibility target.

---

## 6. Frozen legacy restart fixture policy

`tests/cpp_tests/fixtures/legacy_thermo_3d8n_damage.dat` is a **permanent
regression resource**:

* generated from the real legacy thermo-mechanical C++ class before its deletion;
* Hexahedra3D8N;
* local thermal damage with non-trivial committed damage;
* production binary Serializer format;
* must **not** be regenerated using SMA;
* validates the specific historical -> current migration path.

It does not imply arbitrary binary compatibility across compilers/platforms/future
Kratos versions.

---

## 7. Automatic restart Properties rebinding

Phase-6C.1 contract:

* **Persistent** serialized state: damage/history, `EquivalentPlasticStrain`
  threshold, flow-rule state, yield/hardening numerical state, nonlocal state per
  existing semantics.
* **Transient** state: `HardeningLaw::mpProperties` (deliberately not serialized).
* At every relevant damage entry point (`CalculateMaterialResponse*`,
  `FinalizeMaterialResponse*`, `CalculateValue`) the Dam law automatically rebinds
  the hardening law from `ConstitutiveLaw::Parameters::GetMaterialProperties()`.

No external/manual `SetProperties`/`ReinitializeMaterialProperties` is required.
Multi-Properties isolation and restart-continuation tests validate this.

---

## 8. `Has()==false` generic-dispatch contract (permanent maintenance rule)

For parameter-aware outputs (specialized Vector/Matrix outputs,
`LOCAL_EQUIVALENT_STRAIN`), the Dam laws intentionally keep
`Has(variable) == false`, so the generic SMA/Dam element integration-point path
reconstructs current kinematics and calls
`CalculateValue(ConstitutiveLaw::Parameters&, ...)`. Setting `Has == true` would
select direct `GetValue` semantics and bypass the required Parameters-based
computation. New Dam constitutive outputs must follow this convention.

---

## 9. Specialized-output semantics

* `THERMAL_STRAIN = epsilon_th`.
* Linear: `THERMAL_STRESS = C*epsilon_th`, `MECHANICAL_STRESS = C*epsilon_total`.
* Damage: `THERMAL_STRESS = (1-d) C epsilon_th`,
  `MECHANICAL_STRESS = (1-d) C epsilon_total`, with the **same** current damage
  factor for both components.
* `TOTAL_STRESS = MECHANICAL_STRESS - THERMAL_STRESS`.
* Matrix outputs derive from vector outputs.
* Dimensions: 3D vector 6 / tensor 3x3; 2D vector 3 / tensor 2x2.
* Shear convention: strain tensor `xy = gamma_xy/2`, stress tensor `xy = tau_xy`.

---

## 10. Intentional thermal-strain behavior correction

The historical legacy thermo element returned the **total strain** for
`THERMAL_STRAIN_VECTOR`. The permanent contract is
`THERMAL_STRAIN_VECTOR/TENSOR = actual thermal strain`. This is an intentional bug
fix, not an accidental migration difference. The old behavior is not supported.

---

## 11. Nodal-stress smoothing architecture

Dam smoothing scheme workflow:

1. reset nodal accumulators;
2. generic element finalization;
3. `DamNodalCauchyStressExtrapolationProcess` (sole owner);
4. normalize nodal values.

Element-owned `SaveGPStress` / `ExtrapolateGPStress` no longer exist. The process
uses the existing Poro extrapolation utilities where applicable. T3/Q4/T4/H8
behavior is verified; unsupported higher-order extrapolation geometries retain the
established skip behavior; supported geometry with an invalid GP count is an error.

---

## 12. Nonlocal orchestration architecture

```
Dam scheme
  -> DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain
  -> generic SMA Element::CalculateOnIntegrationPoints
  -> Dam law parameter-aware CalculateValue(LOCAL_EQUIVALENT_STRAIN)
  -> Poro averaging
  -> NONLOCAL assigned/read
  -> damage response
```

* No legacy element nonlinear hooks remain.
* LOCAL is current-kinematics transient state.
* NONLOCAL is the averaging-driving value.
* Committed damage/history is separate.
* Both Initialize and Finalize nonlinear-iteration stages remain relevant because
  LOCAL changes after displacement update.

`USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN` remains a genuine runtime activation
flag: Dam solvers set it from `nonlocal_damage`, and
`DamNonlocalDamageUtilities` reads it. There is no alternate legacy ownership or
fallback path anymore. The historical name is retained to avoid unnecessary churn.

---

## 13. Coupled workflows (P / U-P / T-U-P)

* **P-only**: independent of solid mechanics; uses `WaveEquationElement`, pressure
  DOFs and `DamPScheme`.
* **U-P**: solid domain = SMA `SmallDisplacement`; acoustic domain = Dam acoustic
  formulation; coupling via `UPCondition` (nodal U/P, generic condition APIs);
  scheme `DamUPScheme`. No solid concrete-type dependency.
* **T-U-P**: sequential — thermal solve -> nodal TEMPERATURE -> U-P
  mechanical/acoustic stage; thermal behavior comes from Dam constitutive laws on
  SMA solid elements.

---

## 14. Remaining unrelated technical debt

Known issues explicitly left out of scope (not blockers):

* pre-existing 2D `UPCondition` uninitialized-entry quirk for particular
  pure-direction coupling configurations;
* plane-strain thermal semantics (`(1+nu)` factor) intentionally unchanged;
* historical naming of `USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN`.

---

## 15. Permanent test architecture

The permanent Dam test suite covers (at least one meaningful oracle per group):

1. mechanical/registration; 2. thermal linear/nodal; 3. local damage;
4. nonlocal damage; 5. specialized outputs; 6. smoothing;
7. nonlocal orchestration; 8. mass/M1; 9. restart/R4/rebinding;
10. construction/selfweight; 11. P/U-P/T-U-P; 12. `.mdpa` compatibility.
