from math import isclose
import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.linear_strain_energy_response_function import LinearStrainEnergyResponseFunction

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class TestLinearStrainEnergyResponseFunction(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.optimization_problem = OptimizationProblem()
        cls.model_part = cls.model.CreateModelPart("Structure")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

       # create the primal analysis execution policy wrapper
        with kratos_unittest.WorkFolderScope("linear_strain_energy_test", __file__):
            # creating the execution policy wrapper
            execution_policy_wrapper_settings = Kratos.Parameters("""{
                "name"    : "primal",
                "module": "KratosMultiphysics.OptimizationApplication.execution_policies",
                "type": "stepping_analysis_execution_policy",
                "settings": {
                    "model_part_names" : ["Structure"],
                    "analysis_module"  : "KratosMultiphysics.StructuralMechanicsApplication",
                    "analysis_type"    : "StructuralMechanicsAnalysis",
                    "analysis_settings": {
                        "@include_json": "primal_parameters.json"
                    }
                },
                "pre_operations"           : [],
                "post_operations"          : [],
                "log_in_file"              : false,
                "log_file_name"            : "structure.log"
            }""")
            cls.execution_policy_decorator = ExecutionPolicyDecorator(cls.model, execution_policy_wrapper_settings, cls.optimization_problem)
            cls.optimization_problem.AddComponent(cls.execution_policy_decorator)

            Kratos.ModelPartIO("Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)

            # creating the response function wrapper
            response_function_settings = Kratos.Parameters("""{
                "evaluated_model_part_names": ["Structure"],
                "primal_analysis_name"      : "primal",
                "perturbation_size"         : 1e-8
            }""")
            cls.response_function: LinearStrainEnergyResponseFunction = LinearStrainEnergyResponseFunction("strain_energy", cls.model, response_function_settings, cls.optimization_problem)
            cls.optimization_problem.AddComponent(cls.response_function)

            cls.execution_policy_decorator.Initialize()
            cls.response_function.Initialize()

            # now replace the properties
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(cls.model["Structure.structure"], cls.model_part.Elements, False)

            cls.execution_policy_decorator.Execute()
            cls.ref_value = cls.response_function.CalculateValue()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("linear_strain_energy_test", __file__):
            DeleteFileIfExisting("Structure.time")

    def _CheckSensitivity(self, response_function, entities, sensitivity_method, update_method, container_expression_data, delta, rel_tol, abs_tol):
        list_of_sensitivities = []
        for entity in entities:
            adjoint_sensitivity = sensitivity_method(entity)
            update_method(entity, delta)
            self.execution_policy_decorator.Execute()
            value = response_function.CalculateValue()
            fd_sensitivity = (value - self.ref_value)/delta
            update_method(entity, -delta)
            self.assertTrue(
                isclose(adjoint_sensitivity, fd_sensitivity, rel_tol=rel_tol, abs_tol=abs_tol),
                msg = f"{adjoint_sensitivity} != {fd_sensitivity} with rel_to = {rel_tol} and abs_tol = {abs_tol}. Difference = {abs(adjoint_sensitivity - fd_sensitivity)}")
            list_of_sensitivities.append(adjoint_sensitivity)

        list_of_sensitivities = numpy.array(list_of_sensitivities)
        self.assertTrue(all(numpy.isclose(list_of_sensitivities, container_expression_data, rtol=rel_tol, atol=abs_tol)))

    def _UpdateProperties(self, variable, entity, delta):
        entity.Properties[variable] += delta

    def _UpdateNodalPositions(self, direction, entity, delta):
        if direction == 0:
            entity.X += delta
            entity.X0 += delta
        if direction == 1:
            entity.Y += delta
            entity.Y0 += delta
        if direction == 2:
            entity.Z += delta
            entity.Z0 += delta

    def test_CalculateValue(self):
        self.assertAlmostEqual(self.ref_value, 71515947.17480606, 6)

    def test_CalculateYoungModulusSensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(self.model_part)])
        self.response_function.CalculateGradient({Kratos.YOUNG_MODULUS: sensitivity})

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.YOUNG_MODULUS_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.YOUNG_MODULUS, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate(),
            1e-7,
            1e-6,
            1e-5)

    def test_CalculatePoissonRatioSensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.ElementExpression(self.model_part)])
        self.response_function.CalculateGradient({Kratos.POISSON_RATIO: sensitivity})

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.Properties[KratosOA.POISSON_RATIO_SENSITIVITY],
            lambda x, y: self._UpdateProperties(Kratos.POISSON_RATIO, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate(),
            1e-7,
            1e-4,
            1e-5)

    def test_CalculateShapeSensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.NodalExpression(self.model_part)])
        self.response_function.CalculateGradient({KratosOA.SHAPE: sensitivity})

        # calculate nodal shape sensitivities
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_X),
            lambda x, y: self._UpdateNodalPositions(0, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 0],
            1e-6,
            1e-4,
            1e-5)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Y),
            lambda x, y: self._UpdateNodalPositions(1, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 1],
            1e-6,
            1e-4,
            1e-5)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Z),
            lambda x, y: self._UpdateNodalPositions(2, x, y),
            sensitivity.GetContainerExpressions()[0].Evaluate()[:, 2],
            1e-6,
            1e-4,
            1e-5)

if __name__ == "__main__":
    kratos_unittest.main()