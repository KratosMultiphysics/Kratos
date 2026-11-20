"""A real isogeometric analysis, so the NURBS gather has control-point data.

`nurbs_sampling` can already sample an IgaApplication geometry exactly and
gather control-point fields through the geometry's own basis, but nothing
here ever ran an IGA ANALYSIS: the gather was only ever fed synthetic
values. This case runs one and hands back the solved model part.

The case is the Scordelis-Lo roof, a shell benchmark that IgaApplication
ships. It is chosen over the application's `single_patch_test` deliberately:
that one is structurally singular (a direct solver reports a zero column,
and the iterative solver it is configured with returns NaN displacements
without complaining), which its own test never notices because it asserts
nothing numerically.

As with the cylinder and contact cases, the fixture is REUSED IN PLACE and
availability-gated rather than copied.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos

_IGA_TESTS = Path(__file__).resolve().parents[3] / "IgaApplication" / "tests"
_CASE_NAME = "scordelis_roof_test/scordelis_roof_shell_3p"
_PARAMETERS_FILE = _IGA_TESTS / (_CASE_NAME + "_project_parameters.json")


def IsAvailable() -> bool:
    """True when IgaApplication's test data is present."""
    return _PARAMETERS_FILE.is_file()


def _ReadParameters(echo_level: int) -> Kratos.Parameters:
    settings = Kratos.Parameters(_PARAMETERS_FILE.read_text())
    settings["solver_settings"]["echo_level"].SetInt(echo_level)
    if settings["problem_data"].Has("echo_level"):
        settings["problem_data"]["echo_level"].SetInt(echo_level)

    # every file the modelers and the material settings name is relative to
    # the owning tests directory
    if settings.Has("modelers"):
        modelers = settings["modelers"]
        for index in range(modelers.size()):
            parameters = modelers[index]["Parameters"]
            for key in ("geometry_file_name", "physics_file_name", "input_filename"):
                if parameters.Has(key):
                    parameters[key].SetString(str(_IGA_TESTS / parameters[key].GetString()))
    solver = settings["solver_settings"]
    if solver.Has("material_import_settings"):
        materials = solver["material_import_settings"]
        if materials.Has("materials_filename"):
            materials["materials_filename"].SetString(
                str(_IGA_TESTS / materials["materials_filename"].GetString()))
    if settings.Has("output_processes"):
        settings.RemoveValue("output_processes")
    # the reference-result comparisons belong to the owning application's
    # regression and want files relative to its own directory
    if settings.Has("processes"):
        processes = settings["processes"]
        for key in processes.keys():
            kept = Kratos.Parameters("[]")
            for index in range(processes[key].size()):
                entry = processes[key][index]
                name = (entry["python_module"].GetString()
                        if entry.Has("python_module") else "")
                if name in ("from_json_check_result_process",
                            "compare_two_files_check_process"):
                    continue
                kept.Append(entry)
            processes[key] = kept
    return settings


def CreateIgaAnalysis(model: Kratos.Model, echo_level: int = 0):
    """A ready-to-Run() isogeometric shell analysis."""
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
        StructuralMechanicsAnalysis)
    import KratosMultiphysics.IgaApplication  # noqa: F401

    if not IsAvailable():
        raise FileNotFoundError(
            f"IgaApplication test data is not at \"{_PARAMETERS_FILE}\"; this case "
            "reuses it rather than copying the geometry, so the tests that need it "
            "skip when it is absent.")
    return StructuralMechanicsAnalysis(model, _ReadParameters(echo_level))


def GetModelPart(model: Kratos.Model) -> Kratos.ModelPart:
    """The analysis model part, whose NODES are the control points."""
    return model["IgaModelPart"]


def ControlPointDisplacements(model_part: Kratos.ModelPart):
    """(N, 3) DISPLACEMENT at the control points.

    Not a field on the geometry: a control point generally does not lie on
    the surface, which is exactly why the isogeometric gather evaluates the
    NURBS basis instead of reading nodes.
    """
    return numpy.array([
        node.GetSolutionStepValue(Kratos.DISPLACEMENT) for node in model_part.Nodes])


def SolveIga(model: Kratos.Model, echo_level: int = 0):
    """Solves the roof and returns (control_point_displacements, model_part)."""
    analysis = CreateIgaAnalysis(model, echo_level)
    analysis.Run()
    model_part = GetModelPart(model)
    return ControlPointDisplacements(model_part), model_part
