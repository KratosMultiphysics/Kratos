import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication
from fluid_dynamics_analysis import FluidDynamicsAnalysis

try:
    import KratosMultiphysics.MeshMovingApplication
    missing_external_dependencies = True
    missing_application = ''
except ImportError as e:
    missing_external_dependencies = False
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''',
                                    '{0}'.format(e)).group(1)

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities

import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipUnless(missing_external_dependencies,"{} is not available".format(missing_application))
class ALEFluidSolverTest(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    def test_ALEFluidSolver(self):
        work_folder = "test_ale_fluid_solver"
        settings_file_name = "ProjectParameters.json"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteFileIfExisting("ale_fluid_test.time")
            kratos_utilities.DeleteFileIfExisting("test_ale_fluid_solver.post.lst")

    def _runTest(self,settings_file_name):
        model = KM.Model()
        with open(settings_file_name,'r') as settings_file:
            settings = KM.Parameters(settings_file.read())

        # to check the results: add output settings block if needed
        if self.print_output:
            settings["output_processes"].AddValue("gid_output", KM.Parameters(R'''[{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "FluidModelPart.fluid_computational_model_part",
                "output_name"            : "ale_fluid_test",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"       : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"          : "time",
                        "output_control_type" : "step",
                        "output_frequency"    : 1.0,
                        "body_output"         : true,
                        "nodal_results"       : ["VELOCITY","PRESSURE","MESH_DISPLACEMENT","MESH_VELOCITY"]
                    }
                }
            }
            }]'''))

        analysis = FluidDynamicsAnalysis(model,settings)
        analysis.Run()

if __name__ == '__main__':
    UnitTest.main()
