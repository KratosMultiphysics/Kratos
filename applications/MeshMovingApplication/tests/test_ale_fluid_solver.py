import KratosMultiphysics as km
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.MeshMovingApplication
try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError:
    have_external_solvers = False

from fluid_dynamics_analysis import FluidDynamicsAnalysis

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

@UnitTest.skipUnless(have_external_solvers,"Missing required application: ExternalSolversApplication")
class ALEFluidSolverTest(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    def testALEFluidSolver(self):
        work_folder = "test_ale_fluid_solver"
        settings_file_name = "ProjectParameters.json"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteFileIfExisting("ale_fluid_test.time")
            kratos_utilities.DeleteFileIfExisting("test_ale_fluid_solver.post.lst")

    def _runTest(self,settings_file_name):
        model = km.Model()
        with open(settings_file_name,'r') as settings_file:
            settings = km.Parameters(settings_file.read())

        # to check the results: add output settings block if needed
        if self.print_output:
            settings["output_processes"].AddValue("gid_output", km.Parameters(R'''[{
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
