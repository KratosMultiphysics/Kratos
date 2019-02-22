import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as UnitTest

import KratosMultiphysics.kratos_utilities as kratos_utilities
fluid_dynamics_is_available = kratos_utilities.IsApplicationAvailable("FluidDynamicsApplication")
if fluid_dynamics_is_available:
    from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipUnless(fluid_dynamics_is_available,"FluidDynamicsApplication is not available")
class ALEFluidSolverTest(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_vtk_output = False
        self.print_gid_output = False

    def test_ALEFluidSolver(self):
        work_folder = "test_ale_fluid_solver"
        settings_file_name = "ProjectParameters_laplacian_fract_step.json"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteFileIfExisting("test_ale_fluid_solver.post.lst")

    def test_ALEFluidSolverOnSubdomains(self):
        work_folder = "test_ale_fluid_solver"
        settings_file_name = "ProjectParameters_ale_on_subdom_struct_sim_mono.json"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteFileIfExisting("test_ale_fluid_solver.post.lst")

    def _runTest(self,settings_file_name):
        model = KM.Model()
        with open(settings_file_name,'r') as settings_file:
            settings = KM.Parameters(settings_file.read())

        # to check the results: add output settings block if needed
        if self.print_vtk_output:
            settings["output_processes"].AddValue("vtk_output", KM.Parameters(R'''[{
                "python_module" : "vtk_output_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "VtkOutputProcess",
                "help"          : "This process writes postprocessing files for Paraview",
                "Parameters"    : {
                    "model_part_name"                    : "FluidModelPart.domain",
                    "output_control_type"                : "step",
                    "output_frequency"                   : 1,
                    "file_format"                        : "binary",
                    "output_precision"                   : 7,
                    "output_sub_model_parts"             : false,
                    "folder_name"                        : "vtk_output",
                    "save_output_files_in_folder"        : true,
                    "nodal_solution_step_data_variables" : ["VELOCITY","PRESSURE","MESH_DISPLACEMENT","MESH_VELOCITY"],
                    "nodal_data_value_variables"         : [],
                    "element_data_value_variables"       : [],
                    "condition_data_value_variables"     : []
                }
            }]'''))

        if self.print_gid_output:
            settings["output_processes"].AddValue("gid_output", KM.Parameters(R'''[{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "FluidModelPart.domain",
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
