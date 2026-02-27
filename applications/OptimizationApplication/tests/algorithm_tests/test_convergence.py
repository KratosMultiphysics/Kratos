import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
import KratosMultiphysics.OptimizationApplication.convergence_criteria.max_iter_conv_criterion as max_iter_conv_criterion
import KratosMultiphysics.OptimizationApplication.convergence_criteria.l2_conv_criterion as l2_conv_criterion
import KratosMultiphysics.OptimizationApplication.convergence_criteria.combined_conv_criterion as combined_conv_criterion
import KratosMultiphysics.OptimizationApplication.convergence_criteria.avg_abs_improvement_conv_criterion as avg_abs_improvement_conv_criterion
import KratosMultiphysics.OptimizationApplication.convergence_criteria.target_value_conv_criterion as target_value_conv_criterion
import KratosMultiphysics.OptimizationApplication.convergence_criteria.magnitude_reduction_conv_criterion as magnitude_reduction_conv_criterion
import KratosMultiphysics.OptimizationApplication.convergence_criteria.constraint_conv_criterion as constraint_conv_criterion
import KratosMultiphysics.OptimizationApplication.convergence_criteria.patience_conv_criterion as patience_conv_criterion

class TestConvergence(kratos_unittest.TestCase):
    @classmethod
    def setUp(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.CreateElements()

        cls.optimization_problem = OptimizationProblem()
        ComponentDataView("algorithm", cls.optimization_problem).SetDataBuffer(10)

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

    def test_MaxIter(self):
        param = Kratos.Parameters("""{
            "max_iter": 2
        }""")
        convergence_criterium = max_iter_conv_criterion.MaxIterConvCriterion(param, self.optimization_problem)
        convergence_criterium.Initialize()
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        self.assertTrue(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        self.assertTrue(convergence_criterium.IsConverged())

    def test_L2(self):
        param = Kratos.Parameters("""{
            "tolerance": 1e-3
        }""")
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        convergence_criterium = l2_conv_criterion.L2ConvCriterion(param, self.optimization_problem)
        convergence_criterium.Initialize()
        search_direction = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(self.model_part.Elements, Kratos.DENSITY)])
        search_direction.CollectData()
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        search_direction.data[:] /= 100000
        self.assertTrue(convergence_criterium.IsConverged())

    def test_L2_max(self):
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)

        param = Kratos.Parameters("""{
            "max_iter"          : 2
        }""")
        convergence_criterium_max_iter = max_iter_conv_criterion.MaxIterConvCriterion(param, self.optimization_problem)
        param = Kratos.Parameters("""{
            "tolerance"         : 1e-4
        }""")
        convergence_criterium_additional = l2_conv_criterion.L2ConvCriterion(param, self.optimization_problem)
        param = Kratos.Parameters("""{
            "operator": "or"
        }""")
        convergence_criterium = combined_conv_criterion.CombinedConvCriterion(self.model, param, self.optimization_problem)
        convergence_criterium.Add(convergence_criterium_max_iter)
        convergence_criterium.Add(convergence_criterium_additional)
        convergence_criterium.Initialize()

        search_direction = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(self.model_part.Elements, Kratos.DENSITY)])
        search_direction.CollectData()
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        self.assertTrue(convergence_criterium.IsConverged())

    def test_AverAbsDelta(self):
        param = Kratos.Parameters("""{
            "tolerance"         : 2e-3
        }""")
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        convergence_criterium = avg_abs_improvement_conv_criterion.AvgAbsImprovementConvCriterion(param, self.optimization_problem)
        convergence_criterium.Initialize()
        algorithm_data.GetBufferedData()["std_obj_value"] = 5
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.2
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.3
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.2
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.3
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.2222
        self.assertTrue(convergence_criterium.IsConverged())

    def test_AverAbsDelta_max(self):
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        param = Kratos.Parameters("""{
            "max_iter"          : 2
        }""")
        convergence_criterium_max_iter = max_iter_conv_criterion.MaxIterConvCriterion(param, self.optimization_problem)
        param = Kratos.Parameters("""{
            "tolerance"         : 2e-3
        }""")
        convergence_criterium_additional = avg_abs_improvement_conv_criterion.AvgAbsImprovementConvCriterion(param, self.optimization_problem)
        param = Kratos.Parameters("""{
            "operator": "or"
        }""")
        convergence_criterium = combined_conv_criterion.CombinedConvCriterion(self.model, param, self.optimization_problem)
        convergence_criterium.Add(convergence_criterium_max_iter)
        convergence_criterium.Add(convergence_criterium_additional)
        convergence_criterium.Initialize()

        algorithm_data.GetBufferedData()["std_obj_value"] = 5
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.2
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 4.3
        self.assertTrue(convergence_criterium.IsConverged())

    def test_TargetValue(self):
        param = Kratos.Parameters("""{
            "target_value"      : 0.001
        }""")
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        convergence_criterium = target_value_conv_criterion.TargetValueConvCriterion(param, self.optimization_problem)
        convergence_criterium.Initialize()
        algorithm_data.GetBufferedData()["std_obj_value"] = 1
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.1
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.01
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.001
        self.assertTrue(convergence_criterium.IsConverged())

    def test_TargetValue_max(self):
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        param = Kratos.Parameters("""{
            "max_iter"          : 2
        }""")
        convergence_criterium_max_iter = max_iter_conv_criterion.MaxIterConvCriterion(param, self.optimization_problem)
        param = Kratos.Parameters("""{
            "target_value"      : 0.001
        }""")
        convergence_criterium_additional = target_value_conv_criterion.TargetValueConvCriterion(param, self.optimization_problem)
        param = Kratos.Parameters("""{
            "operator": "or"
        }""")
        convergence_criterium = combined_conv_criterion.CombinedConvCriterion(self.model, param, self.optimization_problem)
        convergence_criterium.Add(convergence_criterium_max_iter)
        convergence_criterium.Add(convergence_criterium_additional)
        convergence_criterium.Initialize()

        algorithm_data.GetBufferedData()["std_obj_value"] = 1
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.1
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.01
        self.assertTrue(convergence_criterium.IsConverged())

    def test_TargetMagnitudeReduction(self):
        param = Kratos.Parameters("""{
            "target_scaling_factor"           : 1e-2
        }""")
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        convergence_criterium = magnitude_reduction_conv_criterion.MagnitudeReductionConvCriterion(param, self.optimization_problem)
        convergence_criterium.Initialize()
        algorithm_data.GetBufferedData()["std_obj_value"] = 1
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.1
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.02
        self.assertFalse(convergence_criterium.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.001
        self.assertTrue(convergence_criterium.IsConverged())

    def test_ConstraintConvCriterion(self):
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        param = Kratos.Parameters("""{
            "component_name": "algorithm"
        }""")
        convergence_criterion = constraint_conv_criterion.ConstraintConvCriterion(param, self.optimization_problem)
        convergence_criterion.Initialize()
        algorithm_data.GetBufferedData()["std_value"] = 0.1
        self.assertFalse(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_value"] = -0.1
        self.assertTrue(convergence_criterion.IsConverged())

    def test_CombinedConvCriterion(self):
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        params = Kratos.Parameters("""{
            "operator": "or",
            "list_of_convergence_criteria": [
                {
                    "type": "max_iter_conv_criterion",
                    "settings": {
                        "max_iter": 10
                    }
                },
                {
                    "type": "magnitude_reduction_conv_criterion",
                    "settings": {
                        "target_scaling_factor": 1e-2
                    }
                },
                {
                    "type": "combined_conv_criterion",
                    "settings": {
                        "operator": "and",
                        "list_of_convergence_criteria": [
                            {
                                "type": "target_value_conv_criterion",
                                "settings": {
                                    "target_value"  : 1e-9
                                }
                            },
                            {
                                "type": "l2_conv_criterion",
                                "settings": {
                                   "tolerance" : 1e-4
                                }
                            }
                        ]
                    }
                }
            ]
        }""")
        convergence_criterion = combined_conv_criterion.CombinedConvCriterion(self.model, params, self.optimization_problem)
        convergence_criterion.Initialize()

        search_direction = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(self.model_part.Elements, Kratos.DENSITY)])
        search_direction.CollectData()
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.1

        self.assertFalse(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        algorithm_data.GetBufferedData()["std_obj_value"] = 1e-4
        self.assertTrue(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        search_direction.data[:] /= 1e+6
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        algorithm_data.GetBufferedData()["std_obj_value"] = 0.1
        self.assertFalse(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        search_direction.data[:] /= 1e+6
        algorithm_data.GetBufferedData()["search_direction"] = search_direction
        algorithm_data.GetBufferedData()["std_obj_value"] = 1e-9
        self.assertTrue(convergence_criterion.IsConverged())

    def test_PatienceConvCriterion(self):
        algorithm_data = ComponentDataView("algorithm", self.optimization_problem)
        param = Kratos.Parameters("""{
            "minimum_itr" : 3,
            "patience_itr": 2
        }""")
        convergence_criterion = patience_conv_criterion.PatienceConvCriterion(param, self.optimization_problem)
        convergence_criterion.Initialize()

        algorithm_data.GetBufferedData()["std_obj_value"] = 1e-4
        self.assertFalse(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 1e-4 / 5
        self.assertFalse(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 1e-4 / 10
        self.assertFalse(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 1e-4 / 15
        self.assertFalse(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 1e-4 / 15
        self.assertFalse(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 1e-4 / 150
        self.assertFalse(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 1e-4 / 10
        self.assertFalse(convergence_criterion.IsConverged())
        self.optimization_problem.AdvanceStep()
        algorithm_data.GetBufferedData()["std_obj_value"] = 1e-4 / 10
        self.assertTrue(convergence_criterion.IsConverged())


if __name__ == "__main__":
    kratos_unittest.main()