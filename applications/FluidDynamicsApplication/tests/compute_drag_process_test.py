import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication import compute_drag_process


class ComputeDragProcessTest(KratosUnittest.TestCase):

    def tearDown(self):
        if os.path.exists("main_drag.dat"):
            os.remove("main_drag.dat")

    class DummyDragProcess(compute_drag_process.ComputeDragProcess):
        def _GetFileHeader(self):
            return '# HEADER\n# Time Fx Fy Fz\n'

        def _PrintToScreen(self, result_msg):
            KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedDragProcess", "DUMMY DRAG RESULTS:")
            KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedDragProcess", "Current time: " + result_msg)

        def _GetCorrespondingDragForce(self):
            return [1.0, 2.0, 3.0]

    def testProperFlushing(self):
        model = KratosMultiphysics.Model()
        mpart = model.CreateModelPart("main")

        settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "main",
            "print_drag_to_screen"      : false,
            "write_drag_output_file"    : true,
            "output_file_settings" : {
                "file_extension": "dat",
                "file_name": "main_drag.dat",
                "write_buffer_size": -1
            }
        }""")

        process = self.DummyDragProcess(model, settings)

        process.ExecuteInitialize()

        # After initialization
        with open("main_drag.dat") as f:
            number_of_lines_printed = len(f.readlines())
            self.assertEqual(number_of_lines_printed, 2)

        # First step
        process.ExecuteFinalizeSolutionStep()
        with open("main_drag.dat") as f:
            number_of_lines_printed = len(f.readlines())
            self.assertEqual(number_of_lines_printed, 3)

        # Intermediate step
        process.ExecuteFinalizeSolutionStep()
        with open("main_drag.dat") as f:
            number_of_lines_printed = len(f.readlines())
            self.assertEqual(number_of_lines_printed, 4)

        # End of program
        process.ExecuteFinalize()
        with open("main_drag.dat") as f:
            number_of_lines_printed = len(f.readlines())
            self.assertEqual(number_of_lines_printed, 4)

if __name__ == '__main__':
    KratosUnittest.main()
