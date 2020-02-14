import sys

import KratosMultiphysics as km
import KratosMultiphysics.FluidDynamicsApplication as kfd

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities

class FluidAnalysisTest(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    @UnitTest.skipIf(sys.version_info < (3,0), "this test only runs in Python 3")
    def testFluidDynamicsAnalysis(self):
        work_folder = "CylinderTest"
        settings_file_name = "cylinder_fluid_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteFileIfExisting("cylinder_2d.time")

    def testSteadyAnalysisSmall(self):
        work_folder = "Cavity"
        settings_file_name = "steady_cavity5_fluid_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteFileIfExisting("square5.time")

    def testSteadyCavity(self):
        work_folder = "Cavity"
        settings_file_name = "steady_cavity10_fluid_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteFileIfExisting("square10.time")

    def testSteadyCylinder(self):
        work_folder = "CylinderTest"
        settings_file_name = "steady_cylinder_fluid_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteFileIfExisting("cylinder_2d.time")

    def _runTest(self,settings_file_name):
        model = km.Model()
        with open(settings_file_name,'r') as settings_file:
            settings = km.Parameters(settings_file.read())

        # to check the results: add output settings block if needed
        if self.print_output:
            settings.AddValue("output_processes", km.Parameters(r'''{
                "gid_output" : [{
                    "python_module" : "gid_output_process",
                    "kratos_module" : "KratosMultiphysics",
                    "process_name"  : "GiDOutputProcess",
                    "help"          : "This process writes postprocessing files for GiD",
                    "Parameters"    : {
                        "model_part_name"        : "fluid_computational_model_part",
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
                                "output_frequency"    : 1,
                                "body_output"         : true,
                                "node_output"         : false,
                                "skin_output"         : false,
                                "plane_output"        : [],
                                "nodal_results"       : ["VELOCITY","PRESSURE"],
                                "gauss_point_results" : []
                            },
                            "point_data_configuration"  : []
                        }
                    }
                }]
            }'''))

        analysis = FluidDynamicsAnalysis(model,settings)
        analysis.Run()

if __name__ == '__main__':
    UnitTest.main()

