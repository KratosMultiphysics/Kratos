import os
import sys
from unittest.case import skipIf

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.python_solver import PythonSolver

dependencies_are_available = KratosUtils.CheckIfApplicationsAvailable("FluidDynamicsApplication")
if dependencies_are_available:
    import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

class DummyAnalysis(AnalysisStage):
    def _CreateSolver(self):
        return PythonSolver(self.model, self.project_parameters["solver_settings"])

@KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication")
class TestAnalysisStage(KratosUnittest.TestCase):
    def test_GetMapOfAdditionalHistoricalVariables(self):
        model = KratosMultiphysics.Model()
        settings = KratosMultiphysics.Parameters("""{
            "problem_data":{
                "echo_level": 0,
                "parallel_type" : "OMP"
            },
            "solver_settings":{},
            "processes":{
                "initial_conditions":[
                    {
                        "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
                        "python_module" : "apply_inlet_process",
                        "Parameters" :{
                            "model_part_name" : "ModelPart.inlet"
                        }
                    },
                    {
                        "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
                        "python_module" : "apply_inlet_process",
                        "Parameters" :{
                            "model_part_name" : "ModelPart"
                        }
                    }
                ]
            }
        }""")
        stage = DummyAnalysis(model, settings)
        obtained_map = stage._GetMapOfAdditionalHistoricalVariables()
        expected_map = {"ModelPart": {"VELOCITY"}}
        self.assertEqual(obtained_map, expected_map)
        wrong_map = {"ModelPart": {"VELOCITY", "DISTANCE"}}
        self.assertNotEqual(obtained_map, wrong_map)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
