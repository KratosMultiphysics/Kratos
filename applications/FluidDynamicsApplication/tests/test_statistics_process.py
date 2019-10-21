import KratosMultiphysics as km
import KratosMultiphysics.FluidDynamicsApplication as kfd

from fluid_analysis_without_solution import FluidAnalysisWithoutSolution

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities

class IntegrationPointStatisticsTest(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = True

    def testSteadyAnalysisSmall(self):
        work_folder = "Cavity"
        settings_file_name = "statistics_test_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name)

            kratos_utilities.DeleteFileIfExisting("square5.time")

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
                    "Parameters"    : {
                        "model_part_name"        : "MainModelPart.fluid_computational_model_part",
                        "output_name"            : "cavity",
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

        analysis = FluidAnalysisWithoutSolution(model,settings)
        analysis.Run()

if __name__ == '__main__':
    UnitTest.main()
