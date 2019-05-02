from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as UnitTest

# Other imports
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utilities

import os

# Check other applications dependency
hdf5_is_available = kratos_utilities.CheckIfApplicationsAvailable("HDF5Application")
structural_mechanics_is_available = kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication")

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

    def test_LiftAndMoment(self):
        file_name = "naca0012_small"
        settings_file_name = file_name + "_parameters.json"
        work_folder = "naca0012_small_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteTimeFiles(".")


    def test_Naca0012SmallAdjoint(self):
        if not hdf5_is_available:
            self.skipTest("Missing required application: HDF5Application")
        if not structural_mechanics_is_available:
            self.skipTest("Missing required application: StructuralMechanicsApplication")
        file_name = "naca0012_small_sensitivities"
        settings_file_name_primal = file_name + "_primal_parameters.json"
        settings_file_name_adjoint = file_name + "_adjoint_parameters.json"
        work_folder = "naca0012_small_adjoint_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name_primal)
            self._runTest(settings_file_name_adjoint)

            for file_name in os.listdir(os.getcwd()):
                if file_name.endswith(".h5") or file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

    def test_SmallLiftJumpTest(self):
        file_name = "small_lift_jump"
        settings_file_name = file_name + "_parameters.json"
        work_folder = "naca0012_small_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteTimeFiles(".")

    def _runTest(self,settings_file_name):
        model = KratosMultiphysics.Model()
        with open(settings_file_name,'r') as settings_file:
            settings = KratosMultiphysics.Parameters(settings_file.read())

        if self.print_output:
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
                                "nodal_results"       : ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL","DISTANCE"],
                                "nodal_nonhistorical_results": ["TRAILING_EDGE"],
                                "elemental_conditional_flags_results": ["STRUCTURE"],
                                "gauss_point_results" : ["PRESSURE","VELOCITY","VELOCITY_LOWER","PRESSURE_LOWER","WAKE","ELEMENTAL_DISTANCES","KUTTA"]
                            },
                            "point_data_configuration"  : []
                        }
                    }
                }]
            }'''))

        potential_flow_analysis = PotentialFlowAnalysis(model, settings)
        potential_flow_analysis.Run()

if __name__ == '__main__':
    UnitTest.main()
