from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

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

# Check other applications dependency
hdf5_is_available = kratos_utilities.CheckIfApplicationsAvailable("HDF5Application")
meshing_is_available = kratos_utilities.CheckIfApplicationsAvailable("MeshingApplication")
try:
    import stl
    numpy_stl_is_available = True
except:
    numpy_stl_is_available = False

class WorkFolderScope:
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

    def test_Naca0012SmallAdjoint(self):
        if not hdf5_is_available:
            self.skipTest("Missing required application: HDF5Application")
        file_name = "naca0012_small_sensitivities"
        settings_file_name_primal = file_name + "_primal_parameters.json"
        settings_file_name_adjoint = file_name + "_adjoint_parameters.json"
        settings_file_name_adjoint_analytical = file_name + "_adjoint_analytical_parameters.json"
        work_folder = "naca0012_small_adjoint_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name_primal)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.327805503865, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.MOMENT_COEFFICIENT], -0.105810071870, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], 0.3230253050805644, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.32651526722535246, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.DRAG_COEFFICIENT_FAR_FIELD], 0.0036897206842046205, 0.0, 1e-9)
            self._runTest(settings_file_name_adjoint)
            self._runTest(settings_file_name_adjoint_analytical)

            for file_name in os.listdir(os.getcwd()):
                if file_name.endswith(".h5"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

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

            for file_name in os.listdir():
                if file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)
    def test_EmbeddedCircleNoWake(self):
        if not meshing_is_available:
            self.skipTest("Missing required application: MeshingApplication")
        settings_file_name = "embedded_circle_no_wake_parameters.json"
        work_folder = "embedded_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)

    def test_EmbeddedCircle(self):
        if not hdf5_is_available:
            self.skipTest("Missing required application: HDF5Application")
        if not meshing_is_available:
            self.skipTest("Missing required application: MeshingApplication")
        settings_file_name = "embedded_circle_parameters.json"
        settings_adjoint_file_name = "embedded_circle_adjoint_parameters.json"
        settings_penalty_file_name = "embedded_circle_penalty_parameters.json"
        work_folder = "embedded_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], -0.08769331821378197, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], -0.5405047994795951, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], -0.04198874676923284, 0.0, 1e-9)
            self._runTest(settings_adjoint_file_name)
            self._runTest(settings_penalty_file_name)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.12976919914058177, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], -0.4636936459965071, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.08091879125682809, 0.0, 1e-9)

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
            self._runTest(settings_file_name)
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
            self._runTest(settings_file_name)
            reference_wake_elements_id_list = [1, 2, 3, 4, 5, 6]
            self._validateWakeProcess(reference_wake_elements_id_list, "WAKE")
            reference_kutta_elements_id_list = [13, 14, 15, 17, 18, 19, 20, 21, 22, 23]
            self._validateWakeProcess(reference_kutta_elements_id_list, "KUTTA")

    def test_WakeProcess3DKuttaNodesAboveTheWake(self):
        if not numpy_stl_is_available:
            self.skipTest("Missing required dependency: numpy-stl.")
        # This tests a model with some kutta nodes above the wake
        settings_file_name = "small_3d_parameters.json"
        work_folder = "wake_process_3d_tests/24_elements_kutta_node_above_wake_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)
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
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.7331131286069874, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.DRAG_COEFFICIENT], 0.06480686535448453, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], 0.7228720706323188, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.7287060122732945, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.DRAG_COEFFICIENT_FAR_FIELD], 0.008517301562764179, 0.0, 1e-9)

    def _validateWakeProcess(self,reference_element_id_list, variable_name):
        variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
        solution_element_id_list = []
        for elem in self.main_model_part.Elements:
            if(elem.GetValue(variable)):
                solution_element_id_list.append(elem.Id)
        self._validateIdList(solution_element_id_list, reference_element_id_list)

    def _validateIdList(self, solution_element_id_list, reference_element_id_list):
        if(abs(len(reference_element_id_list) - len(solution_element_id_list)) > 0.1):
            raise Exception('Lists have different lengths', ' reference_element_id_list = ',
                            reference_element_id_list, ' solution_element_id_list = ', solution_element_id_list)
        else:
            for i in range(len(reference_element_id_list)):
                self._check_results(solution_element_id_list[i], reference_element_id_list[i], 0.0, 1e-9)

    def _runTest(self,settings_file_name):
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
                                    "output_frequency"    : 1,
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
        potential_flow_analysis.Run()
        self.main_model_part = model.GetModelPart(settings["solver_settings"]["model_part_name"].GetString())

    def _check_results(self, result, reference, rel_tol, abs_tol):
        isclosethis = t_isclose(result, reference, rel_tol, abs_tol)

        full_msg =  "Failed with following parameters:\n"
        full_msg += str(result) + " != " + str(reference) + ", rel_tol = "
        full_msg += str(rel_tol) + ", abs_tol = " + str(abs_tol)

        self.assertTrue(isclosethis, msg=full_msg)

if __name__ == '__main__':
    UnitTest.main()
