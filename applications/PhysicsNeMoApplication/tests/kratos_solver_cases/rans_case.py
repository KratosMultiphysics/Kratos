"""A turbulent channel flow with a k-epsilon closure.

RANSApplication is compiled in the reference build and nothing here had
driven it. The quantity of interest for a surrogate is the modelled
TURBULENT_VISCOSITY: a learned closure predicts it from the resolved flow
instead of solving the two transport equations for it.

The application's own parameters file is a TEMPLATE - it carries
placeholders such as <STABILIZATION_METHOD> that its test driver substitutes
per test - so this case substitutes them too rather than shipping a second
copy of the file.

As with the other borrowed fixtures, the data is REUSED IN PLACE and
availability-gated.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos

_RANS_TESTS = (Path(__file__).resolve().parents[3]
               / "RANSApplication" / "tests" / "ChannelFlowTest")
_PARAMETERS_FILE = _RANS_TESTS / "channel_flow_mon_ke_parameters.json"

# what the owning application's driver fills in; the short forms only ever
# name reference files, which this case drops anyway
_SUBSTITUTIONS = {
    "<STABILIZATION_METHOD>": "algebraic_flux_corrected",
    "<SHORT_STABILIZATION_METHOD>": "afc",
    "<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>": "velocity_based",
    "<SHORT_WALL_FRICTION_VELOCITY_CALCULATION_METHOD>": "velocity",
    "<FLOW_SOLVER_FORMULATION>": "vms",
    "<LINEAR_SOLVER_TYPE>": "LinearSolversApplication.sparse_lu",
}
_DROPPED_PROCESSES = ("from_json_check_result_process", "compare_two_files_check_process",
                      "json_output_process")


def IsAvailable() -> bool:
    """True when RANSApplication's channel-flow data is present."""
    return _PARAMETERS_FILE.is_file() and (_RANS_TESTS / "channel_flow.mdpa").is_file()


def _ReadParameters(end_time: float, echo_level: int) -> Kratos.Parameters:
    text = _PARAMETERS_FILE.read_text()
    for placeholder, value in _SUBSTITUTIONS.items():
        text = text.replace(placeholder, value)
    settings = Kratos.Parameters(text)
    settings["problem_data"]["end_time"].SetDouble(end_time)
    if settings["problem_data"].Has("echo_level"):
        settings["problem_data"]["echo_level"].SetInt(echo_level)

    solver = settings["solver_settings"]
    if solver.Has("model_import_settings"):
        imports = solver["model_import_settings"]
        if imports.Has("input_filename"):
            name = Path(imports["input_filename"].GetString()).name
            imports["input_filename"].SetString(str(_RANS_TESTS / name))
    if solver.Has("material_import_settings"):
        materials = solver["material_import_settings"]
        if materials.Has("materials_filename"):
            name = Path(materials["materials_filename"].GetString()).name
            materials["materials_filename"].SetString(str(_RANS_TESTS / name))

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


def CreateRansAnalysis(model: Kratos.Model, end_time: float = 1.0, echo_level: int = 0):
    """A ready-to-Run() k-epsilon channel-flow analysis."""
    from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis

    if not IsAvailable():
        raise FileNotFoundError(
            f"RANSApplication channel-flow data is not at \"{_RANS_TESTS}\"; this case "
            "reuses it rather than copying the mesh, so the tests that need it skip "
            "when it is absent.")
    return RANSAnalysis(model, _ReadParameters(end_time, echo_level))


def GetModelPart(model: Kratos.Model) -> Kratos.ModelPart:
    return model["FluidModelPart"]


def ClosureFields(model_part: Kratos.ModelPart):
    """The resolved flow a closure reads, and the viscosity it predicts.

    Returns:
        (features, turbulent_viscosity): an (N, 4) array of VELOCITY_X,
        VELOCITY_Y, turbulent kinetic energy and its dissipation rate, and
        the (N,) TURBULENT_VISCOSITY the k-epsilon model produced.
    """
    import KratosMultiphysics.RANSApplication  # noqa: F401  (registers the variables)

    # the turbulence quantities are registered by the APPLICATION, not by the
    # core, so they are looked up in the global registry rather than reached
    # through the Kratos namespace
    kinetic_energy = Kratos.KratosGlobals.GetVariable("TURBULENT_KINETIC_ENERGY")
    dissipation = Kratos.KratosGlobals.GetVariable(
        "TURBULENT_ENERGY_DISSIPATION_RATE")
    features = numpy.array([
        [node.GetSolutionStepValue(Kratos.VELOCITY_X),
         node.GetSolutionStepValue(Kratos.VELOCITY_Y),
         node.GetSolutionStepValue(kinetic_energy),
         node.GetSolutionStepValue(dissipation)]
        for node in model_part.Nodes])
    viscosity = numpy.array([
        node.GetSolutionStepValue(
            Kratos.KratosGlobals.GetVariable("TURBULENT_VISCOSITY"))
        for node in model_part.Nodes])
    return features, viscosity


def SolveRans(model: Kratos.Model, end_time: float = 1.0):
    """Solves the channel and returns (features, turbulent_viscosity, part)."""
    analysis = CreateRansAnalysis(model, end_time)
    analysis.Run()
    model_part = GetModelPart(model)
    features, viscosity = ClosureFields(model_part)
    return features, viscosity, model_part
