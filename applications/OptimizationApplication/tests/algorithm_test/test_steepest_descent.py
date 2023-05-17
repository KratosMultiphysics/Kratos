from abc import ABC, abstractmethod

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_steepest_descent_new import KratosSteepestDescent
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control import MaterialPropertiesControl

class TestKratosSteepestDescent(kratos_unittest.TestCase, ABC):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.CreateElements()

        ### create response

        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test"
            ]
        }""")
        cls.response_function = MassResponseFunction("mass", cls.model, default_settings)

        ### create control

        parameters = Kratos.Parameters("""{
            "model_part_names"      : ["test"],
            "control_variable_name" : "DENSITY"
        }""")

        cls.properties_control = MaterialPropertiesControl("control1", cls.model, parameters)

        ### create and fill the optimization problem object

        cls.optimization_problem = OptimizationProblem()
        cls.optimization_problem.AddComponent(cls.response_function)
        cls.optimization_problem.AddComponent(cls.properties_control)

        ### optimization algorithm object

        cls.parameters = Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "model_part_name"   : "OptimizationModelPart",
            "objective"         : {
                "response_name": "mass",
                "type"         : "minimization",
                "scaling"      : 1.0
            },
            "controls"          : ["control1"],
            "echo_level"        : 0,
            "settings"          : {
                "gradient_scaling": "inf_norm",
                "echo_level"      : 0,
                "step_size"       : 0.1, 
                "max_iter"        : 2
            }
        }""")

        cls.algorithm = KratosSteepestDescent(cls.model, cls.parameters, cls.optimization_problem)

    @classmethod
    def CreateElements(cls):
        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        properties = cls.model_part.CreateNewProperties(1)
        properties[Kratos.DENSITY] = 2.0
        properties[Kratos.THICKNESS] = 3.0
        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)

        properties = cls.model_part.CreateNewProperties(2)
        properties[Kratos.DENSITY] = 4.0
        properties[Kratos.THICKNESS] = 6.0
        cls.model_part.CreateNewElement("Element2D3N", 2, [4, 1, 3], properties)

    def test_ComputeSearchDirection(self):
        conv = self.algorithm.SolveOptimizationProblem()
        self.assertTrue(conv)
        self.assertEqual(self.algorithm.GetOptimizedObjectiveValue(), 14.249999999999998)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()