from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as UnitTest

# Other imports
from potential_flow_analysis import PotentialFlowAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utilities
import os

class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

class PotentialFlowTestFactory(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    def test_execution(self):

        with WorkFolderScope(self.work_folder):
            self._run_test()

            kratos_utilities.DeleteFileIfExisting("naca0012_Case_5.time")

    def _run_test(self):
        model = KratosMultiphysics.Model()
        with open(self.file_name + "_parameters.json",'r') as parameter_file:
            ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

        if self.print_output:
            ProjectParameters.AddValue("output_processes", KratosMultiphysics.Parameters(r'''{
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
                                "nodal_results"       : ["POSITIVE_FACE_PRESSURE","NEGATIVE_FACE_PRESSURE","TEMPERATURE","UPPER_SURFACE","LOWER_SURFACE","DISTANCE","UPPER_WAKE","LOWER_WAKE","POTENTIAL_JUMP","AIRFOIL","TRAILING_EDGE","KUTTA","DEACTIVATED_WAKE","ZERO_VELOCITY_CONDITION"],
                                "gauss_point_results" : ["PRESSURE","VELOCITY","VELOCITY_LOWER","PRESSURE_LOWER","THICKNESS","TRAILING_EDGE","ELEMENTAL_DISTANCES","KUTTA","ALL_TRAILING_EDGE","ZERO_VELOCITY_CONDITION","DISTANCE","DECOUPLED_TRAILING_EDGE_ELEMENT","TRAILING_EDGE_ELEMENT"]
                            },
                            "point_data_configuration"  : []
                        }
                    }
                }]
            }'''))

        potential_flow_analysis = PotentialFlowAnalysis(model,ProjectParameters)
        potential_flow_analysis.Run()

class Naca0012Test(PotentialFlowTestFactory):
    file_name = "naca0012_Case_5"
    work_folder = "naca0012_tests"

if __name__ == '__main__':
    UnitTest.main()
