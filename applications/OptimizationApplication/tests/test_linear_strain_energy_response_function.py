from math import isclose

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerEnum
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function_wrapper import ResponseFunctionWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function_wrapper import CreateResponseFunctionWrapper

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class TestLinearStrainEnergyResponseFunctionBase(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.optimization_info = OptimizationInfo()
        cls.model_part = cls.model.CreateModelPart("Structure")

       # create the primal analysis execution policy wrapper
        with kratos_unittest.WorkFolderScope("linear_element_test", __file__):
            # creating the execution policy wrapper
            execution_policy_wrapper_settings = Kratos.Parameters("""{
                "name"                     : "primal",
                "execution_policy_settings": {
                    "module"  : "KratosMultiphysics.OptimizationApplication.execution_policies",
                    "type"    : "SteppingAnalysisExecutionPolicy",
                    "settings": {
                        "model_part_names" : ["Structure"],
                        "analysis_settings": {
                            "module": "KratosMultiphysics.StructuralMechanicsApplication",
                            "type": "StructuralMechanicsAnalysis",
                            "settings": {
                                "@include_json": "primal_parameters.json"
                            }
                        }
                    }
                },
                "pre_operations"           : [],
                "post_operations"          : [],
                "log_in_file"              : false,
                "echo_level"               : 0
            }""")
            cls.execution_policy_wrapper = ExecutionPolicyWrapper(cls.model, execution_policy_wrapper_settings)
            cls.optimization_info.AddOptimizationRoutine(cls.execution_policy_wrapper)

            Kratos.ModelPartIO("Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)

            # creating the response function wrapper
            response_function_wrapper_settings = Kratos.Parameters("""{
                "name"     : "strain_energy",
                "module"   : "KratosMultiphysics.OptimizationApplication.responses",
                "type"     : "LinearStrainEnergyResponseFunction",
                "objective": "maximization",
                "scaling"  : 2.0,
                "settings" : {
                    "model_part_name"      : "Structure.structure",
                    "primal_analysis_name" : "primal",
                    "perturbation_size"    : 1e-8
                }
            }""")
            cls.response_function_wrapper = CreateResponseFunctionWrapper(cls.model, response_function_wrapper_settings, cls.optimization_info)
            cls.optimization_info.AddOptimizationRoutine(cls.response_function_wrapper)

            cls.execution_policy_wrapper.Initialize()
            cls.response_function_wrapper.Initialize()

            # now replace the properties
            process_parameters = Kratos.Parameters("""{
                "model_part_name": "Structure.structure",
                "variables_list" : ["DENSITY", "YOUNG_MODULUS", "POISSON_RATIO"],
                "echo_level"     : 0
            }""")
            cls.process: Kratos.Process = KratosOA.EntitySpecificPropertiesProcess(cls.model, process_parameters)

            cls.response_function: ResponseFunctionWrapper = cls.optimization_info.GetOptimizationRoutine("ResponseFunctionWrapper", "strain_energy")
            cls.execution_policy_wrapper.InitializeSolutionStep()
            cls.response_function_wrapper.InitializeSolutionStep()
            cls.process.ExecuteInitializeSolutionStep()
            cls.ref_value = cls.response_function.GetValue()
            cls.standardize_ref_value = cls.response_function.GetStandardizedValue()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("linear_element_test", __file__):
            DeleteFileIfExisting("Structure.time")

    def _CalculateSensitivity(self, sensitivity_variable, sensitivity_container_type: ContainerEnum):
        self.response_function.GetSensitivity(sensitivity_variable, self.model_part, sensitivity_container_type)

    def _CheckSensitivity(self, response_function, entities, sensitivity_method, update_method, delta, rel_tol, abs_tol):
        for entity in entities:
            adjoint_sensitivity = sensitivity_method(entity)
            update_method(entity, delta)
            response_function.InitializeSolutionStep()
            value = response_function.GetValue()
            fd_sensitivity = (value - self.ref_value)/delta
            update_method(entity, -delta)
            self.assertTrue(
                isclose(adjoint_sensitivity, fd_sensitivity, rel_tol=rel_tol, abs_tol=abs_tol),
                msg = f"{adjoint_sensitivity} != {fd_sensitivity} with rel_to = {rel_tol} and abs_tol = {abs_tol}. Difference = {abs(adjoint_sensitivity - fd_sensitivity)}")

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
        self.assertAlmostEqual(self.standardize_ref_value, -2 * 71515947.17480606, 6)

    def test_CalculateYoungModulusSensitivity(self):
        self._CalculateSensitivity(KratosOA.YOUNG_MODULUS_SENSITIVITY, ContainerEnum.ELEMENT_PROPERTIES)

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.GetValue(KratosOA.YOUNG_MODULUS_SENSITIVITY),
            lambda x, y: self._UpdateProperties(Kratos.YOUNG_MODULUS, x, y),
            1e-7,
            1e-6,
            1e-5)

    def test_CalculatePoissonRatioSensitivity(self):
        self._CalculateSensitivity(KratosOA.POISSON_RATIO_SENSITIVITY, ContainerEnum.ELEMENT_PROPERTIES)

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.GetValue(KratosOA.POISSON_RATIO_SENSITIVITY),
            lambda x, y: self._UpdateProperties(Kratos.POISSON_RATIO, x, y),
            1e-7,
            1e-4,
            1e-5)

    def test_CalculateDensitySensitivity(self):
        self._CalculateSensitivity(KratosOA.DENSITY_SENSITIVITY, ContainerEnum.ELEMENT_PROPERTIES)

        # calculate element density sensitivity
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Elements,
            lambda x: x.GetValue(KratosOA.DENSITY_SENSITIVITY),
            lambda x, y: self._UpdateProperties(Kratos.DENSITY, x, y),
            1.0,
            1e-6,
            1e-5)

    def test_CalculateShapeSensitivity(self):
        self._CalculateSensitivity(Kratos.SHAPE_SENSITIVITY, ContainerEnum.NODES)
        # calculate nodal shape sensitivities
        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_X),
            lambda x, y: self._UpdateNodalPositions(0, x, y),
            1e-6,
            1e-4,
            1e-5)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Y),
            lambda x, y: self._UpdateNodalPositions(1, x, y),
            1e-6,
            1e-4,
            1e-5)

        self._CheckSensitivity(
            self.response_function,
            self.model_part.Nodes,
            lambda x: x.GetValue(Kratos.SHAPE_SENSITIVITY_Z),
            lambda x, y: self._UpdateNodalPositions(2, x, y),
            1e-6,
            1e-4,
            1e-5)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()