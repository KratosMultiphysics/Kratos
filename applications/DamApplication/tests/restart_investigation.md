# Dam restart investigation

## Core mechanisms and other applications

Core `RestartUtility` saves and loads a serialized ModelPart with `FileSerializer`,
sets `IS_RESTARTED` and `LOAD_RESTART` after loading, and supports time- and
step-labelled files. `save_restart_process.py` wraps that utility as an output
process. `PythonSolver._ImportModelPart` delegates `input_type: rest` to the
utility. Core `AnalysisStage` reads the restored `ProcessInfo[TIME]` when
restarted. Standard GiD output inspects `IS_RESTARTED` when setting up its list
files, although output naming and single-file behavior still require verification
for the intended continuation workflow.

StructuralMechanics delegates import to `PythonSolver` and skips its normal
post-mdpa model preparation and buffer filling on restart. SMA's base solid
serializes its integration method and law vector and skips material
reinitialization under `IS_RESTARTED`. GeoMechanics also uses the base importer
and guards model preparation. Poromechanics delegates import to the base class
and has a restart guard around parts of preparation; its other ProcessInfo and
buffer handling should be checked individually rather than assumed identical to
StructuralMechanics. These standard mechanisms do not require a Dam-specific
serializer.

## Dam integration limitations

Both Dam mechanical and thermomechanical `ImportModelPart` methods currently
accept only mdpa and otherwise raise `Other input options are not yet
implemented.` An actual DamAnalysis probe with `input_type: rest` reaches that
exception before loading. `SaveRestart` is a placeholder in both solvers; the
standard save process may be used as an output process, but that alone does not
provide a loading path or a complete restart lifecycle.

A small change to import alone is unsafe. `_ExecuteAfterReading` creates
computing submodelparts, assigns constitutive laws and sets buffer size, and
also initializes Python solver attributes such as the computing model-part
names. A restart path must restore the latter without rebuilding or reassigning
loaded model data. DamAnalysis initializes process objects and fills its buffer
by resetting TIME and cloning solution steps from configured start time. It also
sets `self.time` from configuration rather than loaded ProcessInfo and can invoke
previous-results mapping. These operations must be guarded or explicitly
validated for restart; the serialized STEP, historical buffer and load states
must remain intact. Process `ExecuteInitialize` and `ExecuteBeforeSolutionLoop`
may themselves change loaded fields, so each enabled process needs a defined
restart contract. The temporal loop uses Dam's custom time advancement rather
than the base AnalysisStage method; STEP and output scheduling need numerical
verification.

DamAnalysis also calls `CleanPreviousFiles` when initializing standard or
legacy GiD output. A restarted run must not erase original outputs: output
naming, scheduling, whether the original files are retained, and behavior for
single- versus multiple-file modes need explicit design and regression tests.
Thermomechanical restart additionally needs restored thermal and mechanical
computing part names, convection-diffusion settings and both strategies, with no
extra strategy initialization calls. Constructor/initialization settings that
write ProcessInfo need assessment against the serialized values. Construction,
nonlocal orchestration, previous-results transfer and output-device processes
are outside the current validated restart scope.

## Numerical evidence

Two focused tests call the standard file-based `RestartUtility` on a 2D Q4
model after gravity. They verify restored TIME = 4, STEP = 4, displacements,
Gauss-point stresses and assembled residual; after strategy initialization they
verify no stress reset, and compare continuation with unchanged and then
increased horizontal body force. A second test uses the Dam thermoelastic law
with a nonuniform nodal temperature and reference-temperature field; both
fields and the mechanical response survive a file restart and an imposed
subsequent temperature change and horizontal load. Vector comparisons require
both a maximum absolute difference <= 1e-10 and, when the reference vector scale
exceeds 1e-10, relative difference <= 1e-10. These are model/strategy tests:
the temperature evolution is prescribed, not solved by a thermal diffusion
strategy, and neither test runs a restarted DamAnalysis.

Existing element-level tests exercise serialization and continuation of
selected Dam local and nonlocal Simo-Ju damage laws. The nonlocal test explicitly
supplies its nonlocal driving quantity, so it does not prove reconstruction of
neighborhoods, orchestration, or a full model's history. Exact continuation for
other constitutive laws and complete coupled models remains unverified.

## Proposed minimal next step

Implement and test the mechanical rest path first, reusing Core
`RestartUtility` and the standard output process. Restore solver-side computing
model-part names without executing mdpa preparation; preserve serialized TIME,
STEP and buffer; avoid restart-time mapping and buffer filling; define
restart-safe process initialization and non-destructive output naming. A real
DamAnalysis gravity test must compare a file restart with an uninterrupted solve
before enabling the feature. Then independently address the coupled thermal
solver with a nonzero evolving thermal field and verify all relevant history.
This work spans import, analysis lifecycle, process and output handling, so the
current evidence does not justify an abbreviated production patch. Prescribed
Gauss-point stress import with a new displacement datum remains separate.

The original `TestTranferUtility2D.gid` input was rerun from a temporary copy;
its GiD result file was byte-identical to the preceding temporary reference.
The source files were not modified. The run's known zero-RHS warnings occur in
the stationary thermal solve, as documented in the initial-state investigation.

The Dam small suite passes 140 tests with these additions. The rebuilt Dam C++
fast suite passes 50 of 51; `UPCouplingPreservesExpectedSign` is the documented
pre-existing failure. Existing uncommitted initial-state investigations were
preserved. No Dam production restart support was added.
