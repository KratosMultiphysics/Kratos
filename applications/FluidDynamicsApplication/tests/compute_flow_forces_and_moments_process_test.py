import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication import compute_flow_forces_and_moments_process


class ComputeFlowForcesProcessTest(KratosUnittest.TestCase):

    def tearDown(self):
        if os.path.exists("main_drag.dat"):
            os.remove("main_drag.dat")

    class DummyFlowForcesAndMomentsProcess(compute_flow_forces_and_moments_process.ComputeFlowForcesAndMomentsProcess):
        def _GetFileHeader(self):
            return '# HEADER\n# Time Fx Fy Fz\n'

        def _PrintToScreen(self, result_msg):
            KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedFlowForcesAndMomentsProcess", "DUMMY DRAG RESULTS:")
            KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedFlowForcesAndMomentsProcess", "Current time: " + result_msg)

        def _GetCorrespondingFlowForcesAndMoments(self):
            flow_force = [1.0, 2.0, 3.0]
            flow_moment = [0.1, 0.2, 0.3]
            return flow_force, flow_moment

    def testProperFlushing(self):
        model = KratosMultiphysics.Model()
        mpart = model.CreateModelPart("main")

        settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "main",
            "print_flow_forces_and_moments_to_screen"      : false,
            "write_flow_forces_output_file"    : true,
            "output_file_settings" : {
                "file_extension": "dat",
                "file_name": "main_flow_force.dat",
                "write_buffer_size": -1
            }
        }""")

        process = self.DummyFlowForcesAndMomentsProcess(model, settings)

        process.ExecuteInitialize()

        # After initialization
        with open("main_flow_force.dat") as f:
            number_of_lines_printed = len(f.readlines())
            self.assertEqual(number_of_lines_printed, 2)

        # First step
        process.ExecuteFinalizeSolutionStep()
        with open("main_flow_force.dat") as f:
            number_of_lines_printed = len(f.readlines())
            self.assertEqual(number_of_lines_printed, 3)

        # Intermediate step
        process.ExecuteFinalizeSolutionStep()
        with open("main_flow_force.dat") as f:
            number_of_lines_printed = len(f.readlines())
            self.assertEqual(number_of_lines_printed, 4)

        # End of program
        process.ExecuteFinalize()
        with open("main_flow_force.dat") as f:
            number_of_lines_printed = len(f.readlines())
            self.assertEqual(number_of_lines_printed, 4)

if __name__ == '__main__':
    KratosUnittest.main()
