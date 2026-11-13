"""A contact patch test: two blocks pressed together.

ContactStructuralMechanicsApplication is compiled in the reference build but
nothing here exercised it. This case drives its smallest ALM frictionless
patch test - two stacked blocks, a prescribed displacement, and the contact
pressure that results - which is the quantity a contact surrogate predicts
and the one a warm start has to get roughly right.

As with the cylinder, the fixture is REUSED IN PLACE rather than copied:
the owning application ships the mesh and parameters, and duplicating them
here would mean maintaining a second copy. The case is therefore
availability-gated.

Measured cost on the reference machine: 8 nodes and 0.07 s.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos

_CONTACT_TESTS = (Path(__file__).resolve().parents[3]
                  / "ContactStructuralMechanicsApplication" / "tests")
_CASE = _CONTACT_TESTS / "ALM_frictionless_contact_test_2D" / "hyper_simple_patch_test"
_PARAMETERS_FILE = Path(str(_CASE) + "_parameters.json")


def IsAvailable() -> bool:
    """True when the owning application's test data is present."""
    return _PARAMETERS_FILE.is_file()


def _ReadParameters(echo_level: int) -> Kratos.Parameters:
    settings = Kratos.Parameters(_PARAMETERS_FILE.read_text())
    solver = settings["solver_settings"]
    solver["echo_level"].SetInt(echo_level)
    # every path in the file is relative to the owning tests directory
    if solver.Has("model_import_settings"):
        imports = solver["model_import_settings"]
        if imports.Has("input_filename"):
            name = imports["input_filename"].GetString()
            imports["input_filename"].SetString(str(_CONTACT_TESTS / name))
    if solver.Has("material_import_settings"):
        materials = solver["material_import_settings"]
        if materials.Has("materials_filename"):
            name = materials["materials_filename"].GetString()
            materials["materials_filename"].SetString(str(_CONTACT_TESTS / name))
    # The mesh is imported by a MODELER, and its path sits two levels down
    # in that modeler's own model_import_settings - not in the solver's. It
    # is also an ARRAY here: this patch test is two blocks, imported as two
    # mdpa files and merged.
    if settings.Has("modelers"):
        modelers = settings["modelers"]
        for index in range(modelers.size()):
            parameters = modelers[index]["Parameters"]
            for holder in (parameters,
                           parameters["model_import_settings"]
                           if parameters.Has("model_import_settings") else None):
                if holder is None or not holder.Has("input_filename"):
                    continue
                entry = holder["input_filename"]
                if entry.IsArray():
                    absolute = Kratos.Parameters("[]")
                    for position in range(entry.size()):
                        value = absolute.AddEmptyValue(str(position)) if False else None
                        absolute.Append(Kratos.Parameters(
                            '"%s"' % str(_CONTACT_TESTS / entry[position].GetString())))
                    holder["input_filename"] = absolute
                else:
                    holder["input_filename"].SetString(
                        str(_CONTACT_TESTS / entry.GetString()))
    if settings.Has("output_processes"):
        settings.RemoveValue("output_processes")
    # The reference-result comparisons belong to the owning application's
    # regression and want files relative to its own directory. They are
    # spread over differently named lists ("json_check_process" here), so
    # the filter goes by the PROCESS, not by the list's name.
    for block in ("_json_output_process", "_output_processes"):
        if settings.Has(block):
            settings.RemoveValue(block)
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
                if name in ("from_json_check_result_process",
                            "compare_two_files_check_process",
                            "json_output_process"):
                    continue
                kept.Append(entry)
            processes[key] = kept
    return settings


def CreateContactAnalysis(model: Kratos.Model, echo_level: int = 0):
    """A ready-to-Run() contact analysis of the patch test."""
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
        StructuralMechanicsAnalysis)
    import KratosMultiphysics.ContactStructuralMechanicsApplication  # noqa: F401

    if not IsAvailable():
        raise FileNotFoundError(
            f"ContactStructuralMechanics test data is not at \"{_PARAMETERS_FILE}\"; "
            "this case reuses it rather than copying the mesh, so the tests that need "
            "it skip when it is absent.")
    return StructuralMechanicsAnalysis(model, _ReadParameters(echo_level))


def GetModelPart(model: Kratos.Model) -> Kratos.ModelPart:
    return model["Structure"]


def ContactPressures(model_part: Kratos.ModelPart):
    """(N,) nodal contact pressure, the field a surrogate predicts."""
    import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

    return numpy.array([
        node.GetSolutionStepValue(CSMA.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE)
        for node in model_part.Nodes])


def Displacements(model_part: Kratos.ModelPart):
    """(N, 3) nodal DISPLACEMENT."""
    return numpy.array([
        node.GetSolutionStepValue(Kratos.DISPLACEMENT) for node in model_part.Nodes])


def SolveContact(model: Kratos.Model, echo_level: int = 0):
    """Solves the patch test and returns (contact_pressures, model_part)."""
    analysis = CreateContactAnalysis(model, echo_level)
    analysis.Run()
    model_part = GetModelPart(model)
    return ContactPressures(model_part), model_part
