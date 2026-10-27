"""Minimal analysis fixture for execution-backend tests.

Not a real solver: "solving" sets PRESSURE = alpha * x on a small generated
model part, so results are exactly predictable from the sample parameters.
Only depends on the Kratos core (and DatasetExportProcess when the optional
"dataset_export" block is present, used by the subprocess-backend template).
"""

import KratosMultiphysics as Kratos


def Create(model, parameters):
    return DummyAnalysis(model, parameters)


class DummyAnalysis:
    def __init__(self, model, parameters):
        self.model = model
        self.parameters = parameters

    def Run(self):
        alpha = self.parameters["dummy_settings"]["alpha"].GetDouble()
        number_of_nodes = self.parameters["dummy_settings"]["number_of_nodes"].GetInt()

        model_part = self.model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        for i in range(number_of_nodes):
            node = model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, alpha * node.X)
        model_part.ProcessInfo[Kratos.STEP] = 1
        model_part.ProcessInfo[Kratos.TIME] = 1.0

        if self.parameters.Has("dataset_export"):
            from KratosMultiphysics.PhysicsNeMoApplication.processes.export.dataset_export_process import DatasetExportProcess
            export = DatasetExportProcess(self.model, self.parameters["dataset_export"])
            export.ExecuteInitialize()
            export.ExecuteFinalizeSolutionStep()
