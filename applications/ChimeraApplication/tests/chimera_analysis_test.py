import time
import KratosMultiphysics as km
import KratosMultiphysics.FluidDynamicsApplication as kfd
import KratosMultiphysics.ChimeraApplication as kchim
try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError:
    have_external_solvers = False

from fluid_chimera_analysis import FluidChimeraAnalysis

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities

@UnitTest.skipUnless(have_external_solvers,"Missing required application: ExternalSolversApplication")
class ChimeraAnalysisTest(UnitTest.TestCase):

    def setUp(self):
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        self.print_output = False


    def test_MonolithicFlowOverCylinder(self):
        start = time.clock()
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test

        work_folder = "chimera_monolithic_simple_test"
        settings_file_name = "test_chimera_monolithic_simple_ProjectParameters.json"
        with UnitTest.WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
        end = time.clock()
        print("Time taken for monolithic chimera simulation",end-start)

    def test_FractionalStepFlowOverCylinder(self):
        start = time.clock()
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        work_folder = "chimera_fractionalstep_simple_test"
        settings_file_name = "test_chimera_fractionalstep_simple_ProjectParameters.json"
        with UnitTest.WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
        end = time.clock()
        print("Time taken for fractional step chimera simulation",end-start)

    def test_MultipleOverlappingPatchMonolithic(self):
        start = time.clock()
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        work_folder = "chimera_multiple_overlapping_patch"
        settings_file_name = "test_chimera_multiple_overlapping_simple_ProjectParameters.json"
        with UnitTest.WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
        end = time.clock()
        print("Time taken for Multiple overlapping chimera simulation",end-start)

    def _run_test(self,settings_file_name):
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
    test_case = ChimeraAnalysisTest()
    test_case.setUp()
    test_case.test_MonolithicFlowOverCylinder()
    test_case.test_FractionalStepFlowOverCylinder()
    test_case.test_MultipleOverlappingPatchMonolithic()
    print("completed all tests")