"""One-dimensional consolidation: a coupled geomechanics transient.

GeoMechanicsApplication is compiled in the reference build and nothing here
had ever driven it. Consolidation is the natural first case: a saturated
column is loaded, the excess pore water pressure decays as water drains,
and the quantity a surrogate predicts - WATER_PRESSURE - evolves over a
sequence of stages.

The application ships this benchmark as ELEVEN staged parameter files, each
continuing the previous one on the same model. That staging is the case's
time axis, so the series here is built by running consecutive stages on one
model rather than by stepping a single analysis.

As with the other borrowed fixtures, the data is REUSED IN PLACE and
availability-gated rather than copied.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos

_GEO_TESTS = (Path(__file__).resolve().parents[3]
              / "GeoMechanicsApplication" / "tests" / "one_dimensional_consolidation")

_DROPPED_PROCESSES = ("from_json_check_result_process", "compare_two_files_check_process",
                      "json_output_process")


def IsAvailable() -> bool:
    """True when GeoMechanicsApplication's consolidation data is present."""
    return (_GEO_TESTS / "ProjectParameters_stage1.json").is_file()


def _ReadParameters(stage: int, echo_level: int) -> Kratos.Parameters:
    path = _GEO_TESTS / f"ProjectParameters_stage{stage}.json"
    if not path.is_file():
        raise FileNotFoundError(f"No consolidation stage {stage} at \"{path}\".")
    settings = Kratos.Parameters(path.read_text())
    solver = settings["solver_settings"]
    if solver.Has("echo_level"):
        solver["echo_level"].SetInt(echo_level)
    if settings["problem_data"].Has("echo_level"):
        settings["problem_data"]["echo_level"].SetInt(echo_level)

    # every path in the staged files is relative to their own directory
    if solver.Has("model_import_settings"):
        imports = solver["model_import_settings"]
        if imports.Has("input_filename"):
            name = Path(imports["input_filename"].GetString()).name
            imports["input_filename"].SetString(str(_GEO_TESTS / name))
    if solver.Has("material_import_settings"):
        materials = solver["material_import_settings"]
        if materials.Has("materials_filename"):
            name = Path(materials["materials_filename"].GetString()).name
            materials["materials_filename"].SetString(str(_GEO_TESTS / name))

    if settings.Has("output_processes"):
        settings.RemoveValue("output_processes")
    if settings.Has("processes"):
        processes = settings["processes"]
        for key in processes.keys():
            if not processes[key].IsArray():
                continue
            kept = Kratos.Parameters("[]")
            for index in range(processes[key].size()):
                entry = processes[key][index]
                name = (entry["python_module"].GetString()
                        if entry.Has("python_module") else "")
                if name in _DROPPED_PROCESSES:
                    continue
                kept.Append(entry)
            processes[key] = kept
    return settings


def CreateConsolidationAnalysis(model: Kratos.Model, stage: int = 1, echo_level: int = 0):
    """A ready-to-Run() analysis of one consolidation stage."""
    from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import (
        GeoMechanicsAnalysis)

    if not IsAvailable():
        raise FileNotFoundError(
            f"GeoMechanics consolidation data is not at \"{_GEO_TESTS}\"; this case "
            "reuses it rather than copying the mesh, so the tests that need it skip "
            "when it is absent.")
    return GeoMechanicsAnalysis(model, _ReadParameters(stage, echo_level))


def GetModelPart(model: Kratos.Model) -> Kratos.ModelPart:
    return model["PorousDomain"]


def WaterPressures(model_part: Kratos.ModelPart):
    """(N,) nodal WATER_PRESSURE - the field consolidation is about."""
    return numpy.array([
        node.GetSolutionStepValue(Kratos.WATER_PRESSURE) for node in model_part.Nodes])


def SolveStages(model: Kratos.Model, stages=(1, 2, 3), echo_level: int = 0):
    """Runs consecutive stages on ONE model, collecting a state per stage.

    The staging is this benchmark's time axis: each stage continues the
    previous one, so the collected pressures are a transient series a
    temporal surrogate can be trained on.

    Returns:
        (len(stages), N) array of nodal water pressures, and the model part.
    """
    states = []
    model_part = None
    for stage in stages:
        analysis = CreateConsolidationAnalysis(model, stage, echo_level)
        analysis.Run()
        model_part = GetModelPart(model)
        states.append(WaterPressures(model_part))
    return numpy.stack(states), model_part
