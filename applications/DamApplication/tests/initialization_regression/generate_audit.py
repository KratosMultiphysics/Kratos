#!/usr/bin/env python3
"""Generate 2025_process_audit.md and 2025_process_audit.csv.

Data compiled from:
- git diff ece5cfed146..359d1ef414 (PR #13472) for applications/DamApplication
- source tracing of every C++ Dam process, Python wrapper and consumer
  (elements / conditions / constitutive laws)
- git log for later fixes (notably #14617)

The CSV mirrors the Markdown table column for column.
"""

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))

COLUMNS = [
    "Process",
    "Python wrapper",
    "Physical role",
    "Variable(s) written",
    "Entity",
    "Historical callback",
    "Current callback",
    "Lifecycle class (L0/L1/L2/L3)",
    "Historical assignment mechanism",
    "Current assignment mechanism",
    "is_fixed/constrained semantics",
    "Called before solver.Initialize? legacy",
    "Called before solver.Initialize? current",
    "Updated each solution step?",
    "Consumer of value",
    "Initialization consumer?",
    "Existing test coverage",
    "Later fixes",
    "Risk",
    "Evidence level (N/P/S)",
    "current_preloop_call_count_before_first_solve",
    "same_time_idempotent",
    "idempotence_reason",
    "state_accumulated",
    "risk_if_non_idempotent",
    "Recommended action",
]

YES = ("YES (equivalent solver-visible state before first solve; the current "
       "implementation additionally executes the callback once more, idempotently, "
       "after solver.Initialize())")
LEGACY = "YES (pre-solver-init via ExecuteInitialize)"
DOUBLE_CALL = "2 (line 235 pre + line 308 post solver.Initialize) — identical TIME=0"

IDEM_OVERWRITE = "overwrites a deterministic nodal value (FastGetSolutionStepValue(var) = value); no accumulation/randomness/entity/table mutation"
IDEM_NONE = "none (deterministic overwrite; same result after one vs two calls)"

# Shared strings
MECH_CONSUMER = (
    "mechanical elements/constitutive law during CalculateLocalSystem/"
    "CalculateMaterialResponse (solve); not during Initialize"
)
THERMAL_CONSUMER = (
    "thermal elements (EulerianConvDiff/ConvDiff2D) in CalculateLocalSystem and "
    "thermomechanical constitutive law (thermal strain) in CalculateMaterialResponse "
    "(solve); not during Initialize"
)
FACE_HEAT_CONSUMER = (
    "thermal face condition (ConvectionDiffusion) reads nodal FACE_HEAT_FLUX during "
    "CalculateLocalSystem (solve); not during Initialize"
)
PRESS_CONSUMER = (
    "structural surface/line load conditions read nodal POSITIVE_FACE_PRESSURE during "
    "CalculateLocalSystem (solve); not during Initialize"
)
CL_CONSUMER = (
    "thermomechanical constitutive law (thermal strain) in CalculateMaterialResponse "
    "(solve); not during Initialize"
)
NODAL_YOUNG_CONSUMER = (
    "*Nodal constitutive laws (LinearElastic3DLawNodal / ThermalLinearElastic3DLawNodal) "
    "in CalculateMaterialResponse (solve); not during Initialize"
)
ADDED_MASS_CONSUMER = (
    "Dam AddedMassCondition reads nodal ADDED_MASS during CalculateLocalSystem "
    "(mass matrix assembly, solve); not during Initialize"
)

# Columns that are identical for all 19 renamed C++ processes
RENAME_COMMON = {
    "Historical callback": "ExecuteInitialize",
    "Current callback": "ExecuteBeforeSolutionLoop",
    "Lifecycle class (L0/L1/L2/L3)": "L0 (equivalent solver-visible initialization state; "
    "callback-count is NOT literally identical — see idempotence columns)",
    "Historical assignment mechanism": "FastGetSolutionStepValue(var) = value (same body, "
    "method renamed only)",
    "Current assignment mechanism": "FastGetSolutionStepValue(var) = value (identical body)",
    "Called before solver.Initialize? legacy": LEGACY,
    "Called before solver.Initialize? current": YES,
    "Initialization consumer?": "No (only during solve)",
    "current_preloop_call_count_before_first_solve": DOUBLE_CALL,
    "same_time_idempotent": "YES (verified by source trace; representative processes "
    "verified by runtime double-call tests)",
    "state_accumulated": "NO",
    "risk_if_non_idempotent": IDEM_NONE,
}


def rename_row(proc, wrapper, role, variables, entity, updated_step, consumer,
                coverage, risk, action, fix_semantics="unchanged (parameter rename only)",
                evidence="S", idempotence_reason=IDEM_OVERWRITE, overrides=None):
    row = {
        "Process": proc,
        "Python wrapper": wrapper,
        "Physical role": role,
        "Variable(s) written": variables,
        "Entity": entity,
        **RENAME_COMMON,
        "is_fixed/constrained semantics": fix_semantics,
        "Updated each solution step?": updated_step,
        "Consumer of value": consumer,
        "Existing test coverage": coverage,
        "Later fixes": "none",
        "Risk": risk,
        "Evidence level (N/P/S)": evidence,
        "idempotence_reason": idempotence_reason,
        "Recommended action": action,
    }
    if overrides:
        row.update(overrides)
    return row


ROWS = []

# ---------------------------------------------------------------- Family T
ROWS.append(rename_row(
    "DamBofangConditionTemperatureProcess",
    "impose_reservoir_temperature_condition_process.py (BOFANGTEMPERATURE)",
    "thermal boundary condition (reservoir temperature profile)",
    "TEMPERATURE", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=THERMAL_CONSUMER,
    coverage="test_bofang_initialization.py (this branch): analytical + lifecycle + "
             "production-wrapper chain; full thermomechanical experiment (A/B2/C)",
    risk="R0",
    action="keep tests; proven numerically identical pre/post #13472 (A vs B2 vs C)",
    fix_semantics="default false; wrapper FORCES true (AddEmptyValue+SetBool) -> "
                  "Fix(TEMPERATURE) on water nodes + Free each step (same pre/post)",
))

ROWS.append(rename_row(
    "DamReservoirConstantTemperatureProcess",
    "impose_reservoir_temperature_condition_process.py (CONSTANTRESERVOIRTEMPERATURE)",
    "thermal boundary condition (constant reservoir temperature)",
    "TEMPERATURE", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=THERMAL_CONSUMER,
    coverage="construction test (impose_reservoir_temperature_condition_process)",
    risk="R0",
    action="none",
    fix_semantics="default false; wrapper FORCES true -> Fix(TEMPERATURE) + Free each step (same pre/post)",
))

ROWS.append(rename_row(
    "DamReservoirMonitoringTemperatureProcess",
    "impose_reservoir_temperature_condition_process.py (MONITORINGRESERVOIRTEMPERATURE)",
    "thermal boundary condition (measured reservoir temperature)",
    "TEMPERATURE", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=THERMAL_CONSUMER,
    coverage="none",
    risk="R1",
    action="optional small process-level test (assigned value + fixity); lifecycle-equivalent by source",
    fix_semantics="default false; wrapper FORCES true -> Fix(TEMPERATURE) + Free each step (same pre/post)",
))

ROWS.append(rename_row(
    "DamFixTemperatureConditionProcess",
    "impose_uniform_temperature_process.py; impose_face_heat_flux_process.py (T uniform)",
    "thermal boundary condition (fixed/uniform temperature)",
    "TEMPERATURE", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=THERMAL_CONSUMER,
    coverage="construction test (impose_uniform_temperature_process)",
    risk="R0",
    action="none",
    fix_semantics="default false; wrapper passes settings value; Fix + Free each step (same pre/post)",
))

ROWS.append(rename_row(
    "DamAzenhaHeatSourceProcess (DamAzenhaHeatFluxProcess)",
    "impose_heat_source_process.py (AzenhaHeatFlux...)",
    "thermal source (hydration heat, Azenha)",
    "HEAT_FLUX, ALPHA_HEAT_SOURCE", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep; tracks ALPHA_HEAT_SOURCE history)",
    consumer="thermal elements (ConvDiff2D) volume source in CalculateLocalSystem (solve)",
    coverage="none",
    risk="R1",
    action="optional focused unit test; lifecycle-equivalent by source (bodies identical)",
))

ROWS.append(rename_row(
    "DamNoorzaiHeatSourceProcess (DamNoorzaiHeatFluxProcess)",
    "impose_heat_source_process.py (NoorzaiHeatFlux...)",
    "thermal source (Noorzai)",
    "HEAT_FLUX", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep; internal call renamed consistently)",
    consumer="thermal elements (ConvDiff2D) volume source in CalculateLocalSystem (solve)",
    coverage="none",
    risk="R1",
    action="optional focused unit test; lifecycle-equivalent by source",
))

ROWS.append(rename_row(
    "DamTSolAirHeatFluxProcess (DamTSolAirHeatFluxProcess)",
    "impose_face_heat_flux_process.py",
    "thermal boundary condition (T-sol-air ambient heat flux)",
    "FACE_HEAT_FLUX", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=FACE_HEAT_CONSUMER,
    coverage="construction test (impose_face_heat_flux_process)",
    risk="R0",
    action="none",
))

ROWS.append(rename_row(
    "DamNodalReferenceTemperatureProcess",
    "impose_nodal_reference_temperature_process.py",
    "thermomechanical reference quantity (nodal reference temperature)",
    "NODAL_REFERENCE_TEMPERATURE", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=CL_CONSUMER,
    coverage="none",
    risk="R1",
    action="optional focused unit test (value present before first solve)",
))

ROWS.append(rename_row(
    "DamGroutingReferenceTemperatureProcess",
    "impose_grouting_reference_temperature_process.py",
    "thermomechanical reference quantity (grouting)",
    "NODAL_REFERENCE_TEMPERATURE", "nodes",
    updated_step="NO (set at init to mInitialValue; updated only at grouting time in "
                 "ExecuteFinalizeSolutionStep)",
    consumer=CL_CONSUMER,
    coverage="none",
    risk="R1",
    action="optional focused unit test; lifecycle-equivalent by source",
))

ROWS.append(rename_row(
    "DamTemperatureByDeviceProcess",
    "impose_temperature_by_device_process.py (temperature_by_device_list)",
    "thermal boundary condition (device-recorded temperature)",
    "TEMPERATURE", "nodes",
    updated_step="YES (only at ExecuteInitializeSolutionStep)",
    consumer=THERMAL_CONSUMER,
    coverage="none",
    risk="R1",
    action="optional focused unit test",
    # this one was NOT renamed (both hist/cur only have ExecuteInitializeSolutionStep)
    overrides={
        "Historical callback": "ExecuteInitializeSolutionStep (no init callback; unchanged)",
        "Current callback": "ExecuteInitializeSolutionStep (no init callback; unchanged)",
        "Called before solver.Initialize? legacy": "N/A (per-step assignment only)",
        "Called before solver.Initialize? current": "N/A (per-step assignment only)",
    },
))

# ---------------------------------------------------------------- Family M
ROWS.append(rename_row(
    "DamNodalYoungModulusProcess",
    "impose_nodal_young_modulus_process.py",
    "material parameter (nodal Young modulus)",
    "NODAL_YOUNG_MODULUS", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=NODAL_YOUNG_CONSUMER,
    coverage="none",
    risk="R1",
    action="optional focused unit test; note Fix() on a non-DOF creates a DOF via "
           "Node::Fix (no equation impact; identical pre/post)",
    fix_semantics="default false; wrapper FORCES true -> Fix(NODAL_YOUNG_MODULUS) "
                  "(Node::Fix creates a DOF; harmless, identical pre/post)",
))

ROWS.append(rename_row(
    "DamInputTableNodalYoungModulusProcess",
    "impose_input_table_nodal_young_modulus_process.py",
    "material parameter (table nodal Young modulus)",
    "NODAL_YOUNG_MODULUS", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep, table by node id)",
    consumer=NODAL_YOUNG_CONSUMER,
    coverage="none",
    risk="R1",
    action="optional focused unit test",
))

ROWS.append(rename_row(
    "DamChemoMechanicalAgingYoungProcess",
    "impose_chemo_mechanical_aging_process.py",
    "material parameter (aging Young modulus)",
    "NODAL_YOUNG_MODULUS", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep; recomputed from time)",
    consumer=NODAL_YOUNG_CONSUMER,
    coverage="none",
    risk="R1",
    action="optional focused unit test; note log(time) -> -inf at t=0 (NaN transient "
           "overwritten each step; identical pre/post)",
))

ROWS.append(rename_row(
    "DamRandomFieldsVariableProcess",
    "impose_2d_random_fields_variable_process.py / impose_3d_random_fields_variable_process.py",
    "random field (material/thermal property)",
    "user-specified (e.g. NODAL_YOUNG_MODULUS, CONDUCTIVITY)", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep, table by node id)",
    consumer="depends on the variable (material or thermal consumers, solve-time)",
    coverage="none",
    risk="R1",
    action="optional focused unit test",
))

# ---------------------------------------------------------------- Family L
ROWS.append(rename_row(
    "DamHydroConditionLoadProcess",
    "impose_water_loads_condition_process.py (HydroLinePressure2D/HydroSurfacePressure3D)",
    "static load (hydrostatic pressure)",
    "POSITIVE_FACE_PRESSURE", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=PRESS_CONSUMER,
    coverage="construction test (impose_water_loads_condition_process)",
    risk="R0",
    action="none",
))

ROWS.append(rename_row(
    "DamUpliftConditionLoadProcess",
    "impose_water_loads_condition_process.py (StraightUplift...)",
    "uplift/hydro load",
    "POSITIVE_FACE_PRESSURE", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=PRESS_CONSUMER,
    coverage="none",
    risk="R1",
    action="optional focused unit test",
))

ROWS.append(rename_row(
    "DamUpliftCircularConditionLoadProcess",
    "impose_water_loads_condition_process.py (CircularUplift...)",
    "uplift/hydro load (circular)",
    "POSITIVE_FACE_PRESSURE", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=PRESS_CONSUMER,
    coverage="none",
    risk="R1",
    action="optional focused unit test",
))

ROWS.append(rename_row(
    "DamWestergaardConditionLoadProcess",
    "impose_water_loads_condition_process.py (HydroDynamic...)",
    "dynamic hydrodynamic pressure (Westergaard)",
    "POSITIVE_FACE_PRESSURE", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=PRESS_CONSUMER,
    coverage="none",
    risk="R1",
    action="optional focused unit test",
))

# ---------------------------------------------------------------- Family D
ROWS.append(rename_row(
    "DamAddedMassConditionProcess",
    "none (instantiated directly; pybind DamAddedMassConditionProcess)",
    "dynamic/added-mass contribution",
    "ADDED_MASS", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep)",
    consumer=ADDED_MASS_CONSUMER,
    coverage="none",
    risk="R1",
    action="treat separately for dynamic runs; lifecycle-equivalent by source "
           "(ADDED_MASS consumed only in CalculateLocalSystem of the added-mass "
           "condition during the solve); add a small process-level test if used",
))

# ---------------------------------------------------------------- Generic table
ROWS.append(rename_row(
    "ApplyComponentTableProcessDam (apply_component_table_process.hpp)",
    "apply_constraint_vector_dam_table_process.py / apply_load_vector_dam_table_process.py "
    "(table branches)",
    "generic table assignment (constraint/load components)",
    "user-specified component (DISPLACEMENT_X, POINT_LOAD_X, ...)", "nodes",
    updated_step="YES (ExecuteInitializeSolutionStep reads table by time)",
    consumer="depends on variable (DOF or load; solve-time)",
    coverage="test_apply_load_vector_dam_processes.py (table path)",
    risk="R1",
    action="keep existing test",
    fix_semantics="default false; Fix only if constrained (same pre/post)",
))

ROWS.append({
    "Process": "ImposeNodalYoungModulusProcess (wrapper)",
    "Python wrapper": "impose_nodal_young_modulus_process.py",
    "Physical role": "material parameter (nodal Young modulus)",
    "Variable(s) written": "NODAL_YOUNG_MODULUS",
    "Entity": "nodes",
    "Historical callback": "ExecuteInitialize",
    "Current callback": "ExecuteBeforeSolutionLoop",
    "Historical assignment mechanism": "DamNodalYoungModulusProcess via wrapper (is_fixed forced true)",
    "Current assignment mechanism": "DamNodalYoungModulusProcess via wrapper (constrained forced true)",
    "is_fixed/constrained semantics": "wrapper forces constrained=true (same pre/post)",
    "Called before solver.Initialize? legacy": LEGACY,
    "Called before solver.Initialize? current": "N/A - wrapper is unimportable",
    "Updated each solution step?": "would be YES (per-step) once importable",
    "Consumer of value": "nodal constitutive laws (solve)",
    "Initialization consumer?": "No",
    "Existing test coverage": "none",
    "Later fixes": "none",
    "Risk": "R3 (discovered during audit; now FIXED on this branch)",
    "Recommended action": "already fixed (class base Process -> KratosMultiphysics.Process); "
                           "covered by the positive test test_nodal_young_modulus",
})

# ---------------------------------------------------------------- Wrapper-only AssignScalarVariableProcess migrations
ROWS.append({
    "Process": "ApplyConstraintVectorDamTableProcess (wrapper)",
    "Python wrapper": "apply_constraint_vector_dam_table_process.py",
    "Physical role": "constraint (prescribed displacement)",
    "Variable(s) written": "DISPLACEMENT_X / _Y / _Z",
    "Entity": "nodes",
    "Historical callback": "ExecuteInitialize",
    "Current callback": "ExecuteBeforeSolutionLoop",
    "Historical assignment mechanism": "ApplyConstantScalarValueProcess(model_part, params), "
                                       "is_fixed passed from JSON (default false)",
    "Current assignment mechanism": "AssignScalarVariableProcess(Model, params), constrained "
                                    "passed from JSON; ApplyComponentTableProcessDam for tables",
    "is_fixed/constrained semantics": "explicitly passed from JSON (intended: true for "
                                      "constraints); identical values pre/post",
    "Called before solver.Initialize? legacy": LEGACY,
    "Called before solver.Initialize? current": YES,
    "Updated each solution step?": "YES (AssignScalarVariableProcess re-applies each step)",
    "Consumer of value": "DISPLACEMENT is the primary DOF; consumed by the solver "
                         "(fixed value on constrained DOFs)",
    "Initialization consumer?": "No (DOF values read at solve; fixity set before first "
                                "SetUpDofSet in both versions)",
    "Existing test coverage": "joint tests (all 8), construction test, bofang case",
    "Later fixes": "none needed (constrained passed explicitly)",
    "Risk": "R0",
    "Recommended action": "none (covered by existing tests)",
})

ROWS.append({
    "Process": "ApplyLoadVectorDamProcess (wrapper)",
    "Python wrapper": "apply_load_vector_dam_process.py",
    "Physical role": "static load (point-load vector)",
    "Variable(s) written": "POINT_LOAD_X / _Y / _Z",
    "Entity": "nodes",
    "Historical callback": "ExecuteInitialize",
    "Current callback": "ExecuteBeforeSolutionLoop",
    "Historical assignment mechanism": "ApplyConstantScalarValueProcess(model_part, params); "
                                       "is_fixed NOT set -> default false -> no fix, value set once",
    "Current assignment mechanism": "AssignScalarVariableProcess(Model, params); post-PR "
                                    "constrained NOT set -> default TRUE -> ApplyFixity on "
                                    "POINT_LOAD_* (no DOF) -> KRATOS_ERROR crash; FIXED by "
                                    "#14617 (explicit constrained=false)",
    "is_fixed/constrained semantics": "default differs! ApplyConstantScalarValueProcess "
                                      "default is_fixed=false; AssignScalarVariableProcess "
                                      "default constrained=true. #14617 sets false explicitly.",
    "Called before solver.Initialize? legacy": LEGACY,
    "Called before solver.Initialize? current": YES,
    "Updated each solution step?": "YES (AssignScalarVariableProcess re-applies each step)",
    "Consumer of value": "structural point-load condition reads POINT_LOAD during "
                         "CalculateLocalSystem (solve)",
    "Initialization consumer?": "No",
    "Existing test coverage": "test_apply_load_vector_dam_processes.py (added by #14617)",
    "Later fixes": "#14617 (2026-07-30): explicit constrained=false for all components",
    "Risk": "R2 historically (real breakage between #13472 and #14617); current R0/R1 (fixed + tested)",
    "Recommended action": "keep regression test; documented in audit",
})

ROWS.append({
    "Process": "ApplyLoadVectorDamTableProcess (wrapper)",
    "Python wrapper": "apply_load_vector_dam_table_process.py",
    "Physical role": "static load (table point-load vector)",
    "Variable(s) written": "POINT_LOAD_X / _Y / _Z",
    "Entity": "nodes",
    "Historical callback": "ExecuteInitialize",
    "Current callback": "ExecuteBeforeSolutionLoop",
    "Historical assignment mechanism": "ApplyConstantScalarValueProcess + "
                                       "ApplyComponentTableProcessDam (tables)",
    "Current assignment mechanism": "AssignScalarVariableProcess (constrained=false after "
                                    "#14617) + ApplyComponentTableProcessDam (tables)",
    "is_fixed/constrained semantics": "same default-difference as ApplyLoadVectorDamProcess; "
                                      "fixed by #14617",
    "Called before solver.Initialize? legacy": LEGACY,
    "Called before solver.Initialize? current": YES,
    "Updated each solution step?": "YES (table read each step)",
    "Consumer of value": "structural point-load condition, CalculateLocalSystem (solve)",
    "Initialization consumer?": "No",
    "Existing test coverage": "test_apply_load_vector_dam_processes.py (added by #14617)",
    "Later fixes": "#14617 (explicit constrained=false)",
    "Risk": "R2 historically; current R0/R1 (fixed + tested)",
    "Recommended action": "keep regression test; documented in audit",
})

ROWS.append({
    "Process": "ImposeThermalParametersScalarValueProcess (wrapper)",
    "Python wrapper": "impose_thermal_parameters_scalar_value_process.py",
    "Physical role": "material/thermal parameter",
    "Variable(s) written": "DENSITY, CONDUCTIVITY, SPECIFIC_HEAT",
    "Entity": "nodes",
    "Historical callback": "ExecuteInitialize",
    "Current callback": "ExecuteBeforeSolutionLoop",
    "Historical assignment mechanism": "ApplyConstantScalarValueProcess(model_part, params) "
                                       "(is_fixed not set -> false)",
    "Current assignment mechanism": "AssignScalarVariableProcess(Model, params) with explicit "
                                    "constrained=false",
    "is_fixed/constrained semantics": "explicit false (not DOFs); equivalent",
    "Called before solver.Initialize? legacy": LEGACY,
    "Called before solver.Initialize? current": YES,
    "Updated each solution step?": "YES",
    "Consumer of value": "thermal elements read nodal CONDUCTIVITY/SPECIFIC_HEAT/DENSITY in "
                         "CalculateLocalSystem (solve)",
    "Initialization consumer?": "No",
    "Existing test coverage": "construction test",
    "Later fixes": "none",
    "Risk": "R0",
    "Recommended action": "none",
})

ROWS.append({
    "Process": "ImposeUniformTemperatureProcess (wrapper)",
    "Python wrapper": "impose_uniform_temperature_process.py",
    "Physical role": "thermal boundary condition (uniform temperature)",
    "Variable(s) written": "TEMPERATURE",
    "Entity": "nodes",
    "Historical callback": "ExecuteInitialize",
    "Current callback": "ExecuteBeforeSolutionLoop",
    "Historical assignment mechanism": "ApplyConstantScalarValueProcess / DamFixTemperatureConditionProcess",
    "Current assignment mechanism": "constrained from settings: true -> "
                                    "DamFixTemperatureConditionProcess; false -> "
                                    "AssignScalarVariableProcess(constrained from settings)",
    "is_fixed/constrained semantics": "passed through; equivalent",
    "Called before solver.Initialize? legacy": LEGACY,
    "Called before solver.Initialize? current": YES,
    "Updated each solution step?": "YES",
    "Consumer of value": "thermal elements + thermo CL (solve)",
    "Initialization consumer?": "No",
    "Existing test coverage": "construction test",
    "Later fixes": "none",
    "Risk": "R0",
    "Recommended action": "none",
})

ROWS.append({
    "Process": "ImposeFaceHeatFluxProcess (wrapper)",
    "Python wrapper": "impose_face_heat_flux_process.py",
    "Physical role": "thermal boundary condition (uniform T + T-sol-air flux)",
    "Variable(s) written": "TEMPERATURE (uniform), FACE_HEAT_FLUX (T-sol-air)",
    "Entity": "nodes",
    "Historical callback": "ExecuteInitialize",
    "Current callback": "ExecuteBeforeSolutionLoop",
    "Historical assignment mechanism": "ApplyConstantScalarValueProcess (is_fixed=false) + "
                                       "DamFixTemperatureConditionProcess + DamTSolAirHeatFluxProcess",
    "Current assignment mechanism": "AssignScalarVariableProcess (constrained=false) + "
                                    "DamFixTemperatureConditionProcess + DamTSolAirHeatFluxProcess",
    "is_fixed/constrained semantics": "explicit false for the uniform-T component; equivalent",
    "Called before solver.Initialize? legacy": LEGACY,
    "Called before solver.Initialize? current": YES,
    "Updated each solution step?": "YES",
    "Consumer of value": "thermal elements + thermal face condition (solve)",
    "Initialization consumer?": "No",
    "Existing test coverage": "construction test",
    "Later fixes": "none",
    "Risk": "R0",
    "Recommended action": "none",
})


def _md_escape(s):
    return s.replace("|", "\\|").replace("\n", " ")


def write_csv(path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for row in ROWS:
            w.writerow({c: row.get(c, "") for c in COLUMNS})


def write_md(path):
    lines = []
    lines.append("# 2025 DamApplication process-migration audit (PR #13472)\n")
    lines.append("Audit of every DamApplication file changed by PR #13472 "
                 "('Dam/assign scalar variable process', merge `359d1ef414`, 2025-06-02, "
                 "pre-PR parent `ece5cfed146`).\n")
    lines.append("Companion CSV: `2025_process_audit.csv` (same columns).\n")
    lines.append("## Method\n")
    lines.append("- Changed-file inventory reconstructed from git (see below).\n")
    lines.append("- Every C++ process: lifecycle callbacks compared historical vs current "
                 "(`git show ece5cfe:<file>` vs current); bodies verified identical except "
                 "the method-name rename (plus `is_fixed`->`constrained` for 7 processes).\n")
    lines.append("- Every Python wrapper: `ExecuteInitialize`->`ExecuteBeforeSolutionLoop` "
                 "forwarding compared historical vs current.\n")
    lines.append("- `dam_analysis.py`: both process-call sites changed "
                 "`ExecuteInitialize`->`ExecuteBeforeSolutionLoop` (selfweight block + main "
                 "block).\n")
    lines.append("- Consumer trace: for each written variable, located the element/condition/"
                 "constitutive law that reads it and verified it reads only during "
                 "`CalculateLocalSystem` / `CalculateMaterialResponse` (solve), never during "
                 "element `Initialize`, `InitializeMaterial`, scheme/strategy `Initialize` or "
                 "model import.\n")
    lines.append("- Later-fix reconstruction: `git log 359d1ef..HEAD` on the affected files "
                 "(found #14617).\n")
    lines.append("\n## Changed-file inventory (PR #13472)\n")
    lines.append("20 C++ processes, 17 Python wrappers, `dam_analysis.py`, 16 test JSON files.\n")
    cpp = [
        "apply_component_table_process.hpp",
        "dam_added_mass_condition_process.hpp",
        "dam_azenha_heat_source_process.hpp",
        "dam_bofang_condition_temperature_process.hpp",
        "dam_chemo_mechanical_aging_young_process.hpp",
        "dam_fix_temperature_condition_process.hpp",
        "dam_grouting_reference_temperature_process.hpp",
        "dam_hydro_condition_load_process.hpp",
        "dam_input_table_nodal_young_modulus_process.hpp",
        "dam_nodal_reference_temperature_process.hpp",
        "dam_nodal_young_modulus_process.hpp",
        "dam_noorzai_heat_source_process.hpp",
        "dam_random_fields_variable_process.hpp",
        "dam_reservoir_constant_temperature_process.hpp",
        "dam_reservoir_monitoring_temperature_process.hpp",
        "dam_t_sol_air_heat_flux_process.hpp",
        "dam_temperature_by_device_process.hpp",
        "dam_uplift_circular_condition_load_process.hpp",
        "dam_uplift_condition_load_process.hpp",
        "dam_westergaard_condition_load_process.hpp",
    ]
    lines.append("**C++ processes (20):** " + ", ".join("`%s`" % c for c in cpp) + "\n")
    py = [
        "apply_constraint_vector_dam_table_process.py",
        "apply_load_vector_dam_process.py",
        "apply_load_vector_dam_table_process.py",
        "dam_analysis.py",
        "impose_2d_random_fields_variable_process.py",
        "impose_3d_random_fields_variable_process.py",
        "impose_chemo_mechanical_aging_process.py",
        "impose_face_heat_flux_process.py",
        "impose_grouting_reference_temperature_process.py",
        "impose_heat_source_process.py",
        "impose_input_table_nodal_young_modulus_process.py",
        "impose_nodal_reference_temperature_process.py",
        "impose_nodal_young_modulus_process.py",
        "impose_reservoir_temperature_condition_process.py",
        "impose_thermal_parameters_scalar_value_process.py",
        "impose_uniform_temperature_process.py",
        "impose_water_loads_condition_process.py",
    ]
    lines.append("**Python wrappers (17):** " + ", ".join("`%s`" % p for p in py) + "\n")
    lines.append("**Tests (16 JSON):** all under `applications/DamApplication/tests/` "
                 "(`construction`, `joint_*` parameters/results).\n")
    lines.append("## Change-type summary (C++ processes)\n")
    lines.append("| Change type | Processes |\n|---|---|\n")
    lines.append("| `is_fixed`->`constrained` AND `ExecuteInitialize`->`ExecuteBeforeSolutionLoop` "
                 "| apply_component_table, dam_bofang, dam_fix_temperature, dam_nodal_young_modulus, "
                 "dam_reservoir_constant, dam_reservoir_monitoring |\n")
    lines.append("| `ExecuteInitialize`->`ExecuteBeforeSolutionLoop` only "
                 "| dam_added_mass, dam_azenha, dam_chemo_mechanical_aging_young, dam_grouting_"
                 "reference_temperature, dam_hydro, dam_input_table_nodal_young_modulus, "
                 "dam_nodal_reference_temperature, dam_random_fields, dam_t_sol_air, "
                 "dam_uplift_circular, dam_uplift, dam_westergaard, dam_noorzai (also renames an "
                 "internal `this->ExecuteInitialize()` call) |\n")
    lines.append("| `is_fixed`->`constrained` only (callback NOT renamed) "
                 "| dam_temperature_by_device (init callback unchanged: only "
                 "`ExecuteInitializeSolutionStep`) |\n")
    lines.append("\n## Lifecycle call-chain analysis\n")
    lines.append("`dam_analysis.py` drives every process via the process factory. Historical vs "
                 "current:\n")
    lines.append("```\nhistorical                 current\n"
                 "solver.AddVariables/Import   (same)\n"
                 "processes constructed        (same)\n"
                 "for p: p.ExecuteInitialize()  for p: p.ExecuteBeforeSolutionLoop()   <- pre solver.Initialize\n"
                 "buffer fill                  (same)\n"
                 "gid_output.Initialize        (same)\n"
                 "solver.Initialize()          solver.Initialize()\n"
                 "for p: p.ExecuteBeforeSolutionLoop()  (historical wrappers: no-op,\n"
                 "                                       only ExecuteInitialize existed)\n"
                 "                                       current: re-applies identical values)\n"
                 "time loop: p.ExecuteInitializeSolutionStep()  (both, per step)\n"
                 "```\n")
    lines.append("Every C++ process (except `dam_temperature_by_device`) overrides the init "
                  "callback, and every wrapper forwards it (historical `ExecuteInitialize`, "
                  "current `ExecuteBeforeSolutionLoop`).\n")
    lines.append("**Refined conclusion (do not claim literal callback-count equivalence):**\n")
    lines.append("- Historical effective lifecycle: `ExecuteInitialize()` -> assignment -> "
                  "`solver.Initialize()` -> standard `ExecuteBeforeSolutionLoop()` (effectively a "
                  "no-op for these wrappers/processes). Assignment executed **once** pre "
                  "solver.Initialize.\n")
    lines.append("- Current effective lifecycle: pre-solver `ExecuteBeforeSolutionLoop()` -> "
                  "assignment -> `solver.Initialize()` -> standard `ExecuteBeforeSolutionLoop()` "
                  "-> **assignment again**. The callback is executed **twice**, both before the "
                  "first solve, at the same TIME (line 235 pre and line 308 post "
                  "solver.Initialize in `dam_analysis.py`).\n")
    lines.append("- Net effect: **equivalent solver-visible initialization state**, with an "
                  "additional idempotent pre-loop execution in the current implementation (see "
                  "idempotence section below). Lifecycle class **L0** refers to this "
                  "state-equivalence, not to equal callback counts.\n")
    lines.append("Because the consumers read the assigned variables only during the solve "
                  "(verified per family below), the lifecycle rename cannot change the numerical "
                  "result (same argument proven experimentally for Bofang: A vs B2 vs C "
                  "bit-identical).\n")
    lines.append("\n## Idempotence of the double pre-loop callback\n")
    lines.append("All affected callbacks were inspected for non-idempotent operations "
                  "(accumulation `+=`, internal counters, random generation, table mutation, "
                  "entity/DOF creation, constitutive internal-variable modification, "
                  "properties/ProcessInfo mutation, irreversible state). **No affected callback "
                  "performs any of these.** Every callback writes the assigned nodal variable "
                  "with `FastGetSolutionStepValue(var) = value`, where `value` is a pure "
                  "deterministic function of (geometry, TIME/DELTA_TIME, input tables, current "
                  "solution-step TEMPERATURE/state).\n")
    lines.append("Details per process family:\n")
    lines.append("- **deterministic overwrite** (geometry/parameters): added_mass, bofang, "
                  "reservoir_constant, reservoir_monitoring, fix_temperature, hydro, uplift, "
                  "uplift_circular, westergaard, nodal_young_modulus, chemo_mechanical_aging, "
                  "grouting_reference_temperature.\n")
    lines.append("- **table read and overwrite** (deterministic for a fixed table): "
                  "apply_component_table (initial value), nodal_reference_temperature, "
                  "input_table_nodal_young_modulus, random_fields (the random field is generated "
                  "**once** in the wrapper `__init__`; the callback only reads the table).\n")
    lines.append("- **recalculate from TIME/DELTA_TIME and overwrite**: noorzai "
                  "(`exp(-alpha*time + 0.5*delta_time)`, same value at the same TIME).\n")
    lines.append("- **reset-then-overwrite** (no accumulation): azenha resets "
                  "`ALPHA_HEAT_SOURCE = mAlphaInitial` on every init call.\n")
    lines.append("- **AssignScalarVariableProcess / ApplyComponentTableProcessDam** (wrappers): "
                  "`SetVariable` overwrite; `ApplyFixity`/`Fix` are idempotent (fixing an already "
                  "fixed DOF is a no-op; `Node::Fix` creates a DOF only on the first call — the "
                  "DOF set is stable across the second call).\n")
    lines.append("- `dam_temperature_by_device` has no init callback (per-step only).\n")
    lines.append("Representative runtime double-call tests (one deterministic thermal process, "
                  "one material-evolution process, one load process, `DamRandomFieldsVariableProcess` "
                  "and `DamAddedMassConditionProcess`) confirm identical assigned values, fixity "
                  "and DOF sets after one vs two calls: "
                  "`tests/test_dam_process_lifecycle.py::test_idempotent_*`.\n")
    lines.append("\n## Assignment-semantics analysis (`ApplyConstantScalarValueProcess` -> "
                  "`AssignScalarVariableProcess`)\n")
    lines.append("| Aspect | ApplyConstantScalarValueProcess (historical) | "
                 "AssignScalarVariableProcess (current) |\n|---|---|---|\n")
    lines.append("| default fixity | `is_fixed: false` | `constrained: true` |\n")
    lines.append("| Fix behaviour | `Node::Fix` (creates a DOF if missing, then fixes) | "
                 "`VariableUtils.ApplyFixity` -> `pGetDof->FixDof` (**throws** if the variable "
                 "has no DOF) |\n")
    lines.append("| Free behaviour | none (applied once) | frees in `ExecuteFinalizeSolutionStep` "
                 "if it was fixed |\n")
    lines.append("| interval | none (always applied) | `IntervalUtility`; default [0,1e30]; "
                 "Dam convention [0,0] = active at t=0 only |\n")
    lines.append("| per-step | no re-apply (applied once at `ExecuteInitialize`) | re-applies at "
                 "every `ExecuteInitializeSolutionStep` inside interval |\n")
    lines.append("| model/modelpart ctor | both | Model only |\n")
    lines.append("| tables | Dam used `ApplyComponentTableProcessDam` for tables | same for Dam; "
                 "AssignScalarVariableProcess also supports csv tables |\n")
    lines.append("**Consequences**\n")
    lines.append("- For Dam processes that pass `constrained` explicitly, behaviour is "
                 "equivalent.\n")
    lines.append("- The **default difference** (`is_fixed=false` vs `constrained=true`) is the "
                 "root cause of the vector-load regression (see below).\n")
    lines.append("- For non-DOF variables with `constrained=true`, historical code silently "
                 "created a DOF and fixed it (`Node::Fix`); current code **throws** "
                 "(`ApplyFixity`). Both are wrong for non-DOF quantities, but only the current "
                 "one surfaces as a hard error.\n")
    lines.append("\n## Known later regression (vector loads)\n")
    lines.append("- **PR #13472** switched `apply_load_vector_dam_process.py` and "
                 "`apply_load_vector_dam_table_process.py` from `ApplyConstantScalarValueProcess` "
                 "(default `is_fixed=false`) to `AssignScalarVariableProcess` without setting "
                 "`constrained` (default `true`).\n")
    lines.append("- Effect: `ApplyFixity(POINT_LOAD_X,...)` runs on a variable that is not a DOF "
                 "in the Dam solvers (`AddDofs` only creates DISPLACEMENT/TEMPERATURE). Empirically "
                 "confirmed on current master: it throws\n")
    lines.append("  `Trying to fix/free dof of variable POINT_LOAD_X but this dof does not exist`.\n")
    lines.append("- **Fix #14617** (2026-07-30, `95e2345ad1a`) added explicit "
                 "`constrained=false` to all components of both load processes and added "
                 "`test_apply_load_vector_dam_processes.py` (values + no DOFs).\n")
    lines.append("- Only #14617 touched those files since #13472: the breakage was live for "
                 "~13 months and is not covered by any older test.\n")
    lines.append("- Remaining analogous uses audited: every current `AssignScalarVariableProcess` "
                 "use in DamApplication either passes `constrained` explicitly or forwards it "
                 "from the settings. No remaining instance of the unguarded default. The "
                 "`apply_constraint_vector_dam_table_process` correctly forwards `constrained` "
                 "(a constraint process).\n")
    lines.append("\n## Newly discovered issue (impose_nodal_young_modulus_process)\n")
    lines.append("- Reproduced on current master: importing "
                 "`KratosMultiphysics.DamApplication.impose_nodal_young_modulus_process` raises\n")
    lines.append("  `NameError: name 'Process' is not defined`.\n")
    lines.append("- Root cause: PR #13472 removed `from KratosMultiphysics import *` from that "
                 "wrapper (and updated `Process.__init__` to `KratosMultiphysics.Process.__init__`) "
                 "but left the class base as bare `Process`, which is now undefined.\n")
    lines.append("- Impact: the wrapper (and therefore the nodal-Young-modulus feature) is "
                 "completely unusable since #13472. No existing test exercises it (confirmed by "
                  "searching all DamApplication tests/ProjectParameters).\n")
    lines.append("- **Resolved on this branch** (see commit '[DamApplication] Fix nodal Young "
                  "modulus process import'): the class base was corrected to "
                  "`KratosMultiphysics.Process` (one line). The regression test was converted "
                  "from the broken-state canary into a positive test "
                  "(`test_dam_process_lifecycle.py::test_nodal_young_modulus`) verifying import, "
                  "production-factory instantiation, Process base class, execution on a minimal "
                  "ModelPart, the expected `NODAL_YOUNG_MODULUS` value, and the documented "
                  "fixity with no unrelated DOF/fixity changes.\n")
    lines.append("\n## Evidence levels\n")
    lines.append("- **N (numerical proof):** full numerical comparison (A vs B2 vs C) "
                  "demonstrates equivalence. Currently only Bofang.\n")
    lines.append("- **P (process-level proof):** source analysis plus a direct process/runtime "
                  "regression establishes identical assigned values, fixity and relevant state, "
                  "and (for the idempotence tests) the same result after one vs two calls.\n")
    lines.append("- **S (source-level equivalence):** source tracing demonstrates equivalence, "
                  "but no dedicated runtime regression exists. **Level S is NOT 'numerically "
                  "proven'.**\n")
    lines.append("Evidence level is recorded per process in the table below.\n")
    lines.append("\n## Architectural note (double pre-loop call)\n")
    lines.append("- The current Dam analysis (`dam_analysis.py`) invokes "
                  "`ExecuteBeforeSolutionLoop()` **twice** around `solver.Initialize()` (line 235 "
                  "pre, line 308 post) for the initialization path.\n")
    lines.append("- All affected processes have been demonstrated idempotent at the identical "
                  "initial state (source trace + representative runtime double-call tests).\n")
    lines.append("- A future cleanup could reconsider this double-call architecture, but it is "
                  "**out of scope** for this branch: `dam_analysis.py` is intentionally left "
                  "unchanged. Any such cleanup must be a separate change with dedicated "
                  "regression tests.\n")
    lines.append("\n## Full audit table\n")
    lines.append("_Risk: R0 negligible / R1 low / R2 medium / R3 high (see decision rules)._\n")
    lines.append("| %s |\n" % " | ".join(COLUMNS))
    lines.append("|%s|\n" % "|".join(["---"] * len(COLUMNS)))
    for row in ROWS:
        lines.append("| %s |\n" % " | ".join(_md_escape(row.get(c, "")) for c in COLUMNS))
    lines.append("\n## Risk summary\n")
    counts = {}
    for row in ROWS:
        r = row["Risk"].split(" ")[0]
        counts[r] = counts.get(r, 0) + 1
    lines.append("| Risk | count |\n|---|---|\n")
    for r in ["R0", "R1", "R2", "R3"]:
        lines.append("| %s | %d |\n" % (r, counts.get(r, 0)))
    lines.append("\n- **R0 (9):** Bofang, reservoir_constant, fix_temperature, TSolAir, hydro, "
                 "apply_constraint_vector (wrapper), impose_thermal_parameters (wrapper), "
                 "impose_uniform_temperature (wrapper), impose_face_heat_flux (wrapper).\n")
    lines.append("- **R1 (15):** reservoir_monitoring, azenha, noorzai, nodal_reference_temperature, "
                 "grouting_reference_temperature, temperature_by_device, nodal_young_modulus, "
                 "input_table_nodal_young_modulus, chemo_mechanical_aging, random_fields, uplift, "
                 "uplift_circular, westergaard, added_mass, apply_component_table. "
                 "Lifecycle-equivalent; several now have process-level runtime regression "
                 "(evidence P: reservoir_monitoring, nodal_reference_temperature, grouting_"
                 "reference_temperature, nodal_young_modulus, input_table_nodal_young_modulus, "
                 "random_fields, added_mass), the rest remain evidence S (source-level only).\n")
    lines.append("- **R2 (2, historical only):** apply_load_vector / apply_load_vector_table — "
                 "real breakage between #13472 and #14617; now fixed and tested.\n")
    lines.append("- **R3 (1):** impose_nodal_young_modulus_process (wrapper) — newly discovered, "
                 "unimportable since #13472; fix reported separately (see dedicated section).\n")
    lines.append("\n## Decisions / recommendations\n")
    lines.append("- No further production-code change is required by this audit for the "
                 "lifecycle/assignment semantics: all are equivalent for the merged behaviour.\n")
    lines.append("- The two real regressions introduced by #13472 are resolved: the vector-load "
                 "fixity bug was fixed by #14617 (with regression tests), and the "
                 "nodal-Young-modulus import failure was fixed on this branch (with a positive "
                 "regression test).\n")
    lines.append("- The double pre-loop callback architecture is intentionally left unchanged "
                 "(see architectural note).\n")
    lines.append("- R1 rows without a dedicated runtime regression (evidence level S) are "
                 "covered by source-level equivalence only; no numerical claim is made for them.\n")
    with open(path, "w") as f:
        f.write("\n".join(lines))


def main():
    lclass_col = "Lifecycle class (L0/L1/L2/L3)"
    evidence_col = "Evidence level (N/P/S)"
    # N = numerical proof (A vs B2 vs C); P = process-level runtime regression;
    # S = source-level equivalence only (NOT numerically proven).
    EVIDENCE = {
        "DamBofangConditionTemperatureProcess": "N",
        "DamReservoirConstantTemperatureProcess": "P",
        "DamReservoirMonitoringTemperatureProcess": "P",
        "DamFixTemperatureConditionProcess": "P",
        "DamTSolAirHeatFluxProcess (DamTSolAirHeatFluxProcess)": "P",
        "DamNodalReferenceTemperatureProcess": "P",
        "DamGroutingReferenceTemperatureProcess": "P",
        "DamNodalYoungModulusProcess": "P",
        "DamInputTableNodalYoungModulusProcess": "P",
        "DamRandomFieldsVariableProcess": "P",
        "DamHydroConditionLoadProcess": "P",
        "DamAddedMassConditionProcess": "P",
        "ApplyConstraintVectorDamTableProcess (wrapper)": "P",
        "ApplyLoadVectorDamProcess (wrapper)": "P",
        "ApplyLoadVectorDamTableProcess (wrapper)": "P",
        "ImposeThermalParametersScalarValueProcess (wrapper)": "P",
        "ImposeUniformTemperatureProcess (wrapper)": "P",
        "ImposeFaceHeatFluxProcess (wrapper)": "P",
        "ImposeNodalYoungModulusProcess (wrapper)": "P",
    }
    for row in ROWS:
        if lclass_col not in row or not row[lclass_col]:
            row[lclass_col] = "L0 (wrapper forwards the renamed callback; value assigned " \
                              "pre solver.Initialize in both versions)"
        if row["Process"] in EVIDENCE:
            row[evidence_col] = EVIDENCE[row["Process"]]
        elif evidence_col not in row or not row[evidence_col]:
            row[evidence_col] = "S"
        for field, default in [
            ("current_preloop_call_count_before_first_solve", DOUBLE_CALL),
            ("same_time_idempotent",
             "YES (SetVariable overwrite + idempotent ApplyFixity; no accumulation)"),
            ("idempotence_reason", IDEM_OVERWRITE),
            ("state_accumulated", "NO"),
            ("risk_if_non_idempotent", IDEM_NONE),
        ]:
            if field not in row or not row[field]:
                row[field] = default
    write_csv(os.path.join(HERE, "2025_process_audit.csv"))
    write_md(os.path.join(HERE, "2025_process_audit.md"))
    print("written", os.path.join(HERE, "2025_process_audit.csv"))
    print("written", os.path.join(HERE, "2025_process_audit.md"))


if __name__ == "__main__":
    main()
