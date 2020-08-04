import KratosMultiphysics as km

import KratosMultiphysics.kratos_utilities as kratos_utilities
hdf5_is_available = kratos_utilities.CheckIfApplicationsAvailable("HDF5Application")

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_analysis import AdjointFluidAnalysis

import KratosMultiphysics.KratosUnittest as UnitTest

import os

@UnitTest.skipUnless(hdf5_is_available, "HDF5Application is not available")
class AdjointFluidTest(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    def testCylinder(self):
        work_folder = "CylinderTest"
        primal_settings_file_name = "cylinder_fluid_parameters.json"
        adjoint_settings_file_name = "cylinder_adjoint_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._run_test(primal_settings_file_name,adjoint_settings_file_name)

            kratos_utilities.DeleteFileIfExisting("cylinder_2d.time")

    def _run_test(self,primal_parameter_file_name,adjoint_parameter_file_name):
        model = km.Model()
        settings = km.Parameters(r'''{}''')

        with open(primal_parameter_file_name,'r') as primal_parameter_file:
            settings.AddValue("primal_settings", km.Parameters(primal_parameter_file.read()))

        with open(adjoint_parameter_file_name,'r') as adjoint_parameter_file:
            settings.AddValue("adjoint_settings", km.Parameters(adjoint_parameter_file.read()))

        # Add hdf5 output to the primal problem
        settings["primal_settings"]["processes"]["auxiliar_process_list"].Append(km.Parameters(r'''{
            "kratos_module" : "KratosMultiphysics.HDF5Application",
            "python_module" : "single_mesh_primal_output_process",
            "Parameters" : {
                "model_part_name" : "MainModelPart",
                "file_settings" : {
                    "file_name" : "primal_output-<time>",
                    "file_access_mode" : "truncate"
                },
                "model_part_output_settings" : {
                    "prefix" : "/ModelData"
                },
                "nodal_solution_step_data_settings" : {
                    "list_of_variables": ["VELOCITY", "ACCELERATION", "PRESSURE"]
                },
                "output_time_settings" : {
                    "step_frequency": 1
                }
            }
        }'''))

        # to check the results: add output settings block if needed
        if self.print_output:
            settings["adjoint_settings"].AddValue("output_processes", km.Parameters(r'''{
                "gid_output" : [{
                    "python_module" : "gid_output_process",
                    "kratos_module" : "KratosMultiphysics",
                    "process_name"  : "GiDOutputProcess",
                    "help"          : "This process writes postprocessing files for GiD",
                    "Parameters"    : {
                        "model_part_name"        : "MainModelPart",
                        "output_name"            : "interface_test",
                        "postprocess_parameters" : {
                            "result_file_configuration" : {
                                "gidpost_flags" : {
                                    "GiDPostMode"           : "GiD_PostBinary",
                                    "WriteDeformedMeshFlag" : "WriteUndeformed",
                                    "WriteConditionsFlag"   : "WriteElementsOnly",
                                    "MultiFileFlag"         : "SingleFile"
                                },
                                "file_label"          : "time",
                                "output_control_type" : "step",
                                "output_interval"     : 1,
                                "body_output"         : true,
                                "node_output"         : false,
                                "skin_output"         : false,
                                "plane_output"        : [],
                                "nodal_results"       : ["VELOCITY","PRESSURE","ADJOINT_FLUID_VECTOR_1","ADJOINT_FLUID_SCALAR_1","SHAPE_SENSITIVITY"],
                                "gauss_point_results" : []
                            },
                            "point_data_configuration"  : []
                        }
                    }
                }]
            }'''))

        primal_analysis = FluidDynamicsAnalysis(model,settings["primal_settings"])
        primal_analysis.Run()
        adjoint_model = km.Model()
        adjoint_analysis = AdjointFluidAnalysis(adjoint_model,settings["adjoint_settings"])
        adjoint_analysis.Run()
        self._remove_h5_files("primal_output")

    def _remove_h5_files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                kratos_utilities.DeleteFileIfExisting(name)

if __name__ == '__main__':
    UnitTest.main()

