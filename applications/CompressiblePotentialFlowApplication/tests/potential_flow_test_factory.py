# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as UnitTest

# Other imports
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.KratosUnittest import isclose as t_isclose

import os

try:
    import stl
    numpy_stl_is_available = True
except:
    numpy_stl_is_available = False

try:
    import KratosMultiphysics.MultilevelMonteCarloApplication
    import xmc
    import exaqute
    import numpy as np
    is_xmc_available = True
except:
    is_xmc_available = False

class WorkFolderScope:
    # TODO use KratosUnittest.WorkFolderScope
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

class PotentialFlowTests(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    @UnitTest.skipIfApplicationsNotAvailable("LinearSolversApplication", "HDF5Application")
    def test_Naca0012SmallAdjoint(self):
        file_name = "naca0012_small_sensitivities"
        settings_file_name_primal = file_name + "_primal_parameters.json"
        settings_file_name_adjoint = file_name + "_adjoint_parameters.json"
        settings_file_name_adjoint_farfield = file_name + "_adjoint_far_field_parameters.json"
        settings_file_name_adjoint_analytical = file_name + "_adjoint_analytical_parameters.json"
        work_folder = "naca0012_small_adjoint_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name_primal)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.327805503865, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.MOMENT_COEFFICIENT], -0.105810071870, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], 0.3230253050805644, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.32651526722535246, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.DRAG_COEFFICIENT_FAR_FIELD], 0.0036897206842046205, 0.0, 1e-9)
            self._check_results(self.main_model_part.GetNode(13).GetValue(CPFApp.POTENTIAL_JUMP), 0.3230253050805686, 0.0, 1e-9)
            self._runTest(settings_file_name_adjoint)
            self._runTest(settings_file_name_adjoint_farfield)
            self._runTest(settings_file_name_adjoint_analytical)

            for file_name in os.listdir(os.getcwd()):
                if file_name.endswith(".h5"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

    @UnitTest.skipIfApplicationsNotAvailable("LinearSolversApplication")
    def test_Naca0012SmallCompressible(self):
        file_name = "naca0012_small_compressible"
        settings_file_name = file_name + "_parameters.json"
        work_folder = "naca0012_small_compressible_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.4968313580730855, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.MOMENT_COEFFICIENT], -0.1631792300021498, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], 0.4876931961465126, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.4953997676243705, 0.0, 1e-9)
            self._check_results(self.main_model_part.GetNode(13).GetValue(CPFApp.POTENTIAL_JUMP), 0.48769319614651147, 0.0, 1e-9)
            self._check_perimeter_computation()

            for file_name in os.listdir():
                if file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

    @UnitTest.skipIfApplicationsNotAvailable("LinearSolversApplication")
    def test_Naca0012SmallTransonic(self):
        file_name = "naca0012_small_transonic"
        settings_file_name = file_name + "_parameters.json"
        work_folder = "naca0012_small_transonic_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.9068468588561012, 0.0, 1e-8)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.MOMENT_COEFFICIENT], -0.3804473187215503, 0.0, 1e-8)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], 0.8890010565994741, 0.0, 1e-8)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.9085576033568474, 0.0, 1e-8)
            self._check_results(self.main_model_part.GetNode(13).GetValue(CPFApp.POTENTIAL_JUMP), 0.88900105659947, 0.0, 1e-9)

        kratos_utilities.DeleteTimeFiles(work_folder)

    @UnitTest.skipIfApplicationsNotAvailable("LinearSolversApplication")
    def test_Naca0012SmallPerturbationCompressible(self):
        file_name = "naca0012_small_perturbation_compressible"
        settings_file_name = file_name + "_parameters.json"
        work_folder = "naca0012_small_perturbation_compressible_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.4968313580730855, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.MOMENT_COEFFICIENT], -0.1631792300021498, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], 0.4876931961465126, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.4953997676243705, 0.0, 1e-9)
            self._check_results(self.main_model_part.GetNode(13).GetValue(CPFApp.POTENTIAL_JUMP), 0.48769319614651147, 0.0, 1e-9)

            for file_name in os.listdir():
                if file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

    @UnitTest.skipIfApplicationsNotAvailable("MeshingApplication")
    def test_EmbeddedCircleNoWake(self):
        settings_file_name = "embedded_circle_no_wake_parameters.json"
        work_folder = "embedded_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)

    @UnitTest.skipIfApplicationsNotAvailable("HDF5Application", "MeshingApplication")
    def test_EmbeddedCircle(self):
        settings_file_name = "embedded_circle_parameters.json"
        settings_adjoint_file_name = "embedded_circle_adjoint_parameters.json"
        settings_adjoint_farfield_file_name = "embedded_circle_adjoint_farfield_parameters.json"
        settings_penalty_file_name = "embedded_circle_penalty_parameters.json"
        settings_comp_penalty_file_name = "embedded_compressible_circle_penalty_parameters.json"
        work_folder = "embedded_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], -1.0503073352419314, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], 0.8267434236568114, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], -1.0964242473327956, 0.0, 1e-9)
            self._runTest(settings_adjoint_file_name)
            self._runTest(settings_adjoint_farfield_file_name)
            self._runTest(settings_penalty_file_name)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.27882962377313825, 0.0, 1e-6)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], -0.38305428616956905, 0.0, 1e-6)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.373749145115962, 0.0, 1e-6)
            self._runTest(settings_comp_penalty_file_name)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.2790356699193499, 0.0, 1e-6)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], -0.3832238355625421, 0.0, 1e-6)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.3739407323566888, 0.0, 1e-6)

            for file_name in os.listdir(os.getcwd()):
                if file_name.endswith(".h5"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

    def test_WakeProcess3DSmall(self):
        if not numpy_stl_is_available:
            self.skipTest("Missing required dependency: numpy-stl.")
        # This tests a simple small 3D model
        settings_file_name = "small_3d_parameters.json"
        work_folder = "wake_process_3d_tests/15_elements_small_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name, initialize_only=True)
            reference_wake_elements_id_list = [2, 4, 9, 13, 15]
            self._validateWakeProcess(reference_wake_elements_id_list, "WAKE")
            reference_kutta_elements_id_list = [1, 10, 14]
            self._validateWakeProcess(reference_kutta_elements_id_list, "KUTTA")

    def test_WakeProcess3DNodesOnWake(self):
        if not numpy_stl_is_available:
            self.skipTest("Missing required dependency: numpy-stl.")
        # This tests a model with nodes laying on the wake
        settings_file_name = "small_3d_parameters.json"
        work_folder = "wake_process_3d_tests/25_elements_nodes_on_wake_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name, initialize_only=True)
            reference_wake_elements_id_list = [13, 14, 15, 16, 17]
            self._validateWakeProcess(reference_wake_elements_id_list, "WAKE")
            reference_kutta_elements_id_list = [18, 19, 20, 21, 22, 23]
            self._validateWakeProcess(reference_kutta_elements_id_list, "KUTTA")

    def test_WakeProcess3DKuttaNodesAboveTheWake(self):
        if not numpy_stl_is_available:
            self.skipTest("Missing required dependency: numpy-stl.")
        # This tests a model with some kutta nodes above the wake
        settings_file_name = "small_3d_parameters.json"
        work_folder = "wake_process_3d_tests/24_elements_kutta_node_above_wake_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name, initialize_only=True)
            reference_wake_elements_id_list = [2, 4, 8, 16, 17, 19, 20]
            self._validateWakeProcess(reference_wake_elements_id_list, "WAKE")
            reference_kutta_elements_id_list = [10, 11, 12, 13, 14, 15, 18, 23, 24]
            self._validateWakeProcess(reference_kutta_elements_id_list, "KUTTA")

    def test_Rhombus3DIncompressible(self):
        if not numpy_stl_is_available:
            self.skipTest("Missing required dependency: numpy-stl.")
        settings_file_name = "rhombus_3d_parameters.json"
        work_folder = "rhombus_3d"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.736320273390592, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[KratosMultiphysics.DRAG_COEFFICIENT], 0.06509077036563415, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], 0.7258065125219674, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.7341491672597356, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.DRAG_COEFFICIENT_FAR_FIELD],0.008730661509157661, 0.0, 1e-9)
            self._check_results(self.main_model_part.GetNode(49).GetValue(CPFApp.POTENTIAL_JUMP), 0.35919110345927313, 0.0, 1e-9)

    @UnitTest.skipIfApplicationsNotAvailable("ShapeOptimizationApplication", "LinearSolversApplication", "MeshMovingApplication")
    def test_ShapeOptimizationLiftConstrainedBodyFitted2D(self):
        work_folder = "body_fitted_opt"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            __import__(work_folder+".run_test")

    @UnitTest.skipIfApplicationsNotAvailable("ShapeOptimizationApplication", "LinearSolversApplication", "MeshMovingApplication", "MultilevelMonteCarloApplication", "MappingApplication")
    def test_StochasticShapeOptimizationLiftConstrainedBodyFitted2D(self):
        if not is_xmc_available:
            self.skipTest("XMC and its dependencies could not be imported. Please check applications/MultilevelMonteCarloApplication/README.md for installation details")

        work_folder = "stochastic_body_fitted_opt"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            __import__(work_folder+".run_test")

    def _validateWakeProcess(self,reference_element_id_list, variable_name):
        variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
        solution_element_id_list = []
        for elem in self.main_model_part.Elements:
            if(elem.GetValue(variable)):
                solution_element_id_list.append(elem.Id)
        self._validateIdList(solution_element_id_list, reference_element_id_list)

    def _validateIdList(self, solution_element_id_list, reference_element_id_list):
        # TODO replace with unittest.assertListEqual or KratosUnitTest.assertVectorAlmostEqual
        if(abs(len(reference_element_id_list) - len(solution_element_id_list)) > 0.1):
            raise Exception('Lists have different lengths', ' reference_element_id_list = ',
                            reference_element_id_list, ' solution_element_id_list = ', solution_element_id_list)
        else:
            for i in range(len(reference_element_id_list)):
                self._check_results(solution_element_id_list[i], reference_element_id_list[i], 0.0, 1e-9)

    def _runTest(self,settings_file_name, initialize_only = False):
        model = KratosMultiphysics.Model()
        with open(settings_file_name,'r') as settings_file:
            settings = KratosMultiphysics.Parameters(settings_file.read())

        if self.print_output:
            if settings_file_name == "naca0012_small_sensitivities_adjoint_parameters.json":
                settings.AddValue("output_processes", KratosMultiphysics.Parameters(r'''{
                    "gid_output" : [{
                        "python_module" : "gid_output_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"  : "GiDOutputProcess",
                        "help"          : "This process writes postprocessing files for GiD",
                        "Parameters"    : {
                            "model_part_name"        : "MainModelPart",
                            "output_name"            : "naca0012_adjoint",
                            "postprocess_parameters" : {
                                "result_file_configuration" : {
                                    "gidpost_flags"       : {
                                        "GiDPostMode"           : "GiD_PostBinary",
                                        "WriteDeformedMeshFlag" : "WriteDeformed",
                                        "WriteConditionsFlag"   : "WriteConditions",
                                        "MultiFileFlag"         : "SingleFile"
                                    },
                                    "file_label"          : "step",
                                    "output_control_type" : "step",
                                    "output_frequency"    : 1,
                                    "body_output"         : true,
                                    "node_output"         : false,
                                    "skin_output"         : false,
                                    "plane_output"        : [],
                                    "nodal_results"       : ["SHAPE_SENSITIVITY","ADJOINT_VELOCITY_POTENTIAL", "ADJOINT_AUXILIARY_VELOCITY_POTENTIAL"],
                                    "nodal_nonhistorical_results": [],
                                    "elemental_conditional_flags_results": [],
                                    "gauss_point_results" : []
                                },
                                "point_data_configuration"  : []
                            }
                        }
                    }]
                }'''))
            else:
                settings.AddValue("output_processes", KratosMultiphysics.Parameters(r'''{
                    "gid_output" : [{
                        "python_module" : "gid_output_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"  : "GiDOutputProcess",
                        "help"          : "This process writes postprocessing files for GiD",
                        "Parameters"    : {
                            "model_part_name"        : "MainModelPart",
                            "output_name"            : "naca0012",
                            "postprocess_parameters" : {
                                "result_file_configuration" : {
                                    "gidpost_flags"       : {
                                        "GiDPostMode"           : "GiD_PostBinary",
                                        "WriteDeformedMeshFlag" : "WriteDeformed",
                                        "WriteConditionsFlag"   : "WriteConditions",
                                        "MultiFileFlag"         : "SingleFile"
                                    },
                                    "file_label"          : "step",
                                    "output_control_type" : "step",
                                    "output_interval"    : 1,
                                    "body_output"         : true,
                                    "node_output"         : false,
                                    "skin_output"         : false,
                                    "plane_output"        : [],
                                    "nodal_results"       : ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"],
                                    "nodal_nonhistorical_results": ["TRAILING_EDGE","WAKE_DISTANCE"],
                                    "elemental_conditional_flags_results": ["STRUCTURE"],
                                    "gauss_point_results" : ["PRESSURE_COEFFICIENT","VELOCITY","WAKE","WAKE_ELEMENTAL_DISTANCES","KUTTA"]
                                },
                                "point_data_configuration"  : []
                            }
                        }
                    }]
                }'''))

        potential_flow_analysis = PotentialFlowAnalysis(model, settings)
        potential_flow_analysis.Initialize()
        if not initialize_only:
            potential_flow_analysis.RunSolutionLoop()
            potential_flow_analysis.Finalize()
        self.main_model_part = model.GetModelPart(settings["solver_settings"]["model_part_name"].GetString())

    def _check_results(self, result, reference, rel_tol, abs_tol):
        # TODO add message directly to t_isclose
        isclosethis = t_isclose(result, reference, rel_tol, abs_tol)

        full_msg =  "Failed with following parameters:\n"
        full_msg += str(result) + " != " + str(reference) + ", rel_tol = "
        full_msg += str(rel_tol) + ", abs_tol = " + str(abs_tol)

        self.assertTrue(isclosethis, msg=full_msg)

    def _check_perimeter_computation(self):
        body_model_part = self.main_model_part.GetSubModelPart("Body2D_Body")
        perimeter1 = sum([cond.GetGeometry().Area() for cond in body_model_part.Conditions])
        perimeter2 = CPFApp.PotentialFlowUtilities.CalculateArea(body_model_part.Conditions)
        self._check_results(perimeter1, perimeter2, 1e-9, 1e-9)

class NightlyPotentialFlowTests(UnitTest.TestCase):

    @UnitTest.skipIfApplicationsNotAvailable("ShapeOptimizationApplication", "LinearSolversApplication", "MeshMovingApplication", "MultilevelMonteCarloApplication", "MappingApplication")
    def test_StochasticShapeOptimizationLiftConstrainedBodyFitted2D(self):
        if not is_xmc_available:
            self.skipTest("XMC and its dependencies could not be imported. Please check applications/MultilevelMonteCarloApplication/README.md for installation details")

        work_folder = "stochastic_body_fitted_opt"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            __import__(work_folder+".run_mlmc_test")

if __name__ == '__main__':
    UnitTest.main()
