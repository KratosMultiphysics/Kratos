"""Runner of the real-solver thermal template case.

Reads case parameters (conductivity/heat flux/divisions) and the dataset
export configuration from ProjectParameters.json, runs an actual
ConvectionDiffusionApplication stationary solve, and exports the solved
TEMPERATURE field with DatasetExportProcess — the results contract of the
active-learning SubprocessBackend.

The thermal_case helper module is made importable through PYTHONPATH by the
caller (the tests add tests/kratos_solver_cases to it).
"""

import KratosMultiphysics as Kratos

import thermal_case

if __name__ == "__main__":
    with open("ProjectParameters.json", "r") as f:
        parameters = Kratos.Parameters(f.read())
    case = parameters["case_settings"]

    model = Kratos.Model()
    analysis = thermal_case.CreateThermalAnalysis(
        model,
        conductivity=case["conductivity"].GetDouble(),
        heat_flux=case["heat_flux"].GetDouble(),
        divisions=case["divisions"].GetInt())
    analysis.Run()

    from KratosMultiphysics.PhysicsNeMoApplication.processes.export.dataset_export_process import DatasetExportProcess
    model_part = model["ThermalModelPart"]
    model_part.ProcessInfo[Kratos.STEP] = 1
    export = DatasetExportProcess(model, parameters["dataset_export"])
    export.ExecuteInitialize()
    export.ExecuteFinalizeSolutionStep()
