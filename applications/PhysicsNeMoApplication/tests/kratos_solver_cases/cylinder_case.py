"""Transient flow past a cylinder: the vortex-shedding case.

The application's own fluid case is a steady lid-driven cavity, which is
the wrong shape for a transient surrogate. This one drives the canonical
cylinder instead, reusing the mesh and parameters that already ship with
FluidDynamicsApplication rather than copying a 216 kB mesh into this
application. The fixture is therefore AVAILABILITY-GATED: if that
application's test data is not present, the case reports so and the tests
that need it skip.

Measured cost on the reference machine: 1110 nodes, 2069 elements, and
0.1 s of wall clock to t = 0.2 - cheap enough for a unit test, which is why
this is a real cylinder rather than a stand-in.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos

_CYLINDER_DIRECTORY = (Path(__file__).resolve().parents[3]
                       / "FluidDynamicsApplication" / "tests" / "CylinderTest")
_PARAMETERS_FILE = _CYLINDER_DIRECTORY / "cylinder_fluid_parameters.json"
_MESH_FILE = _CYLINDER_DIRECTORY / "cylinder_2d.mdpa"

# processes that only compare against reference files: they belong to the
# owning application's own regression, not to a surrogate fixture
_DROPPED_PROCESSES = ("compare_two_files_check_process", "point_output_process")


def IsAvailable() -> bool:
    """True when FluidDynamicsApplication's cylinder test data is present."""
    return _PARAMETERS_FILE.is_file() and _MESH_FILE.is_file()


def _ReadParameters(end_time: float, time_step: float, echo_level: int) -> Kratos.Parameters:
    settings = Kratos.Parameters(_PARAMETERS_FILE.read_text())
    settings["problem_data"]["end_time"].SetDouble(end_time)
    if settings["problem_data"].Has("echo_level"):
        settings["problem_data"]["echo_level"].SetInt(echo_level)
    settings["solver_settings"]["echo_level"].SetInt(echo_level)
    # the mesh path in the file is relative to its own directory
    settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(
        str(_MESH_FILE.with_suffix("")))
    if settings["solver_settings"].Has("time_stepping"):
        stepping = settings["solver_settings"]["time_stepping"]
        if stepping.Has("time_step"):
            stepping["time_step"].SetDouble(time_step)
    # the materials file is named relative to the parameters' own directory
    if settings["solver_settings"].Has("material_import_settings"):
        materials = settings["solver_settings"]["material_import_settings"]
        if materials.Has("materials_filename"):
            name = Path(materials["materials_filename"].GetString()).name
            materials["materials_filename"].SetString(str(_CYLINDER_DIRECTORY / name))
    if settings.Has("output_processes"):
        settings.RemoveValue("output_processes")

    # keep the boundary conditions - dropping the whole processes block
    # leaves the flow undriven, and a case that never moves looks like it ran
    if settings.Has("processes"):
        processes = settings["processes"]
        for key in processes.keys():
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


def CreateCylinderAnalysis(model: Kratos.Model, end_time: float = 0.2,
                           time_step: float = 0.01, echo_level: int = 0):
    """A ready-to-run transient analysis of flow past the cylinder."""
    from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import (
        FluidDynamicsAnalysis)

    if not IsAvailable():
        raise FileNotFoundError(
            f"FluidDynamicsApplication's cylinder test data is not at "
            f"\"{_CYLINDER_DIRECTORY}\"; this case reuses it rather than copying the "
            "mesh, so the tests that need it skip when it is absent.")
    return FluidDynamicsAnalysis(
        model, _ReadParameters(end_time, time_step, echo_level))


def GetModelPart(model: Kratos.Model) -> Kratos.ModelPart:
    """The cylinder analysis's main model part."""
    return model["MainModelPart"]


def CollectVelocities(model_part: Kratos.ModelPart):
    """(N, 2) planar VELOCITY, the shape the transient harness collects."""
    return numpy.array([
        [node.GetSolutionStepValue(Kratos.VELOCITY_X),
         node.GetSolutionStepValue(Kratos.VELOCITY_Y)]
        for node in model_part.Nodes])
