import os
import KratosMultiphysics as km

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.RANSApplication.periodic_fluid_dynamics_analysis import PeriodicFluidDynamicsAnalysis

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities

class EvmKEpsilonTest(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    def testOneElementKEpsilonSteady(self):
        work_folder = "OneElement"
        settings_file_name = "one_element_steady_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name)
            kratos_utilities.DeleteTimeFiles(".")

    def testBackwardFacingStepKEpsilonTransient(self):
        work_folder = "BackwardFacingStepTest"
        settings_file_name = "backward_facing_step_k_epsilon_transient_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name)
            kratos_utilities.DeleteTimeFiles(".")

    @UnitTest.skipUnless(km.IsDistributedRun(), "Running without MPI")
    def testChannelFlowKEpsilonTransientMPI(self):
        work_folder = "ChannelFlowTest"
        settings_file_name = "channel_flow_k_epsilon_transient_mpi_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name)

            for file_name in os.listdir():
                if file_name.startswith("channel_flow_") and file_name.endswith(".mdpa"):
                    kratos_utilities.DeleteFileIfExisting(file_name)
            kratos_utilities.DeleteTimeFiles(".")

    def testChannelFlowKEpsilonSteady(self):
        work_folder = "ChannelFlowTest"
        settings_file_name = "channel_flow_k_epsilon_steady_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name)
            kratos_utilities.DeleteTimeFiles(".")

    @UnitTest.skipUnless(km.IsDistributedRun(), "Running without MPI")
    def testChannelFlowKEpsilonSteadyMPI(self):
        work_folder = "ChannelFlowTest"
        settings_file_name = "channel_flow_k_epsilon_steady_mpi_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name)

            for file_name in os.listdir():
                if file_name.startswith("channel_flow_") and file_name.endswith(".mdpa"):
                    kratos_utilities.DeleteFileIfExisting(file_name)
            kratos_utilities.DeleteTimeFiles(".")

    def testChannelFlowKEpsilonSteadyPeriodic(self):
        work_folder = "ChannelFlowTest"
        settings_file_name = "channel_flow_k_epsilon_steady_periodic_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name, True)
            kratos_utilities.DeleteTimeFiles(".")

    def _runTest(self,settings_file_name, is_periodic = False):
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

        if is_periodic:
            analysis = PeriodicFluidDynamicsAnalysis(model,settings)
        else:
            analysis = FluidDynamicsAnalysis(model,settings)

        analysis.Run()

if __name__ == '__main__':
    UnitTest.main()

