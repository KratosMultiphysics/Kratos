"""A NACA 0012 airfoil in compressible potential flow.

CompressiblePotentialFlowApplication is compiled in the reference build and
`adjoint_bridge` already dispatches its response functions, but nothing here
ever ran one: the aerodynamic path was a dispatch table entry with no case
behind it. This drives the application's own small NACA case, which solves
in hundredths of a second on 22 nodes.

As with the cylinder, contact and IGA cases, the fixture is REUSED IN PLACE
and availability-gated rather than copied.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos

_AERO_TESTS = (Path(__file__).resolve().parents[3]
               / "CompressiblePotentialFlowApplication" / "tests")
_CASE_DIRECTORY = _AERO_TESTS / "naca0012_small_compressible_test"
_PARAMETERS_FILE = _CASE_DIRECTORY / "naca0012_small_compressible_parameters.json"
_MESH_FILE = _AERO_TESTS / "naca0012_small_mdpa" / "naca0012_small.mdpa"

_DROPPED_PROCESSES = ("from_json_check_result_process", "compare_two_files_check_process",
                      "json_output_process")


def IsAvailable() -> bool:
    """True when CompressiblePotentialFlowApplication's test data is present."""
    return _PARAMETERS_FILE.is_file() and _MESH_FILE.is_file()


def _ReadParameters(mach_infinity: float, angle_of_attack, echo_level: int) -> Kratos.Parameters:
    settings = Kratos.Parameters(_PARAMETERS_FILE.read_text())
    settings["solver_settings"]["echo_level"].SetInt(echo_level)
    if settings["problem_data"].Has("echo_level"):
        settings["problem_data"]["echo_level"].SetInt(echo_level)
    settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(
        str(_MESH_FILE.with_suffix("")))

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
                # The flight condition is the case parameter a surrogate
                # learns over. THIS case parameterizes by Mach number only:
                # its far-field process carries no angle_of_attack at all,
                # so setting one would silently do nothing.
                if entry.Has("Parameters"):
                    parameters = entry["Parameters"]
                    if parameters.Has("mach_infinity"):
                        parameters["mach_infinity"].SetDouble(mach_infinity)
                    if angle_of_attack is not None:
                        if not parameters.Has("angle_of_attack"):
                            raise ValueError(
                                "This NACA case has no \"angle_of_attack\" to set; it "
                                "is parameterized by \"mach_infinity\".")
                        parameters["angle_of_attack"].SetDouble(angle_of_attack)
                kept.Append(entry)
            processes[key] = kept
    return settings


def CreateAirfoilAnalysis(model: Kratos.Model, mach_infinity: float = 0.6,
                          angle_of_attack=None, echo_level: int = 0):
    """A ready-to-Run() potential-flow analysis at one flight condition."""
    from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import (
        PotentialFlowAnalysis)

    if not IsAvailable():
        raise FileNotFoundError(
            f"CompressiblePotentialFlow test data is not at \"{_PARAMETERS_FILE}\"; this "
            "case reuses it rather than copying the mesh, so the tests that need it "
            "skip when it is absent.")
    return PotentialFlowAnalysis(
        model, _ReadParameters(mach_infinity, angle_of_attack, echo_level))


def GetModelPart(model: Kratos.Model) -> Kratos.ModelPart:
    return model["FluidModelPart"]


def GetBodyModelPart(model: Kratos.Model) -> Kratos.ModelPart:
    """The airfoil skin - where a surrogate's point cloud lives."""
    return model["FluidModelPart.Body2D_Body"]


def VelocityPotentials(model_part: Kratos.ModelPart):
    """(N,) nodal VELOCITY_POTENTIAL, the solved unknown."""
    import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFA

    return numpy.array([
        node.GetSolutionStepValue(CPFA.VELOCITY_POTENTIAL) for node in model_part.Nodes])


def SolveAirfoil(model: Kratos.Model, mach_infinity: float = 0.6):
    """Solves one flight condition and returns (potentials, model_part)."""
    analysis = CreateAirfoilAnalysis(model, mach_infinity)
    analysis.Run()
    model_part = GetModelPart(model)
    return VelocityPotentials(model_part), model_part
