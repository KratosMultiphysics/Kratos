#    |  /           |
#    ' /   __| _` | __|  _ \   __|
#    . \  |   (   | |   (   |\__ `
#   _|\_\_|  \__,_|\__|\___/ ____/
#                   Multi-Physics
#
#  License:		 BSD License
#					 Kratos default license: kratos/license.txt
#
#  Main authors:    Suneth Warnakulasuriya           
#

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import RetrieveObject

class SteppingAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(model, parameters)

        default_settings = Kratos.Parameters("""{
            "model_part_names" : [],
            "analysis_settings": {
                "module"  : "",
                "type"    : "",
                "settings": {}
            }
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_parts = []
        self.analysis = RetrieveObject(self.model, parameters["analysis_settings"], AnalysisStage)

    def Initialize(self, _: dict):
        self.analysis.Initialize()

        # initialize model parts
        self.model_parts = [self.model[model_part_name] for model_part_name in self.parameters["model_part_names"].GetStringArray()]

    def Execute(self, _: dict):
        time_before_analysis = []
        step_before_analysis = []
        delta_time_before_analysis = []

        for model_part in self.model_parts:
            time_before_analysis.append(model_part.ProcessInfo[Kratos.TIME])
            step_before_analysis.append(model_part.ProcessInfo[Kratos.STEP])
            delta_time_before_analysis.append(model_part.ProcessInfo[Kratos.DELTA_TIME])

        # Reset step/time iterators such that they match the optimization iteration after calling CalculateValue (which internally calls CloneTimeStep)
        for index, model_part in enumerate(self.model_parts):
            model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_analysis[index] - 1)
            model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_analysis[index] - 1)
            model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0)

        self.analysis.time = self.analysis._GetSolver().AdvanceInTime(self.analysis.time)
        self.analysis.InitializeSolutionStep()
        self.analysis._GetSolver().Predict()
        self.analysis._GetSolver().SolveSolutionStep()

        self.analysis.FinalizeSolutionStep()
        self.analysis.OutputSolutionStep()

        # Clear results or modifications on model parts
        for index, model_part in enumerate(self.model_parts):
            model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_analysis[index])
            model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_analysis[index])
            model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, delta_time_before_analysis[index])



