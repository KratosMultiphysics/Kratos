import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as Fluid

import KratosMultiphysics.KratosUnittest as UnitTest

class FlowsMeasuringUtilityTest(UnitTest.TestCase):

    def setUp(self):
        self.current_model = Kratos.Model()
        self.fluid_model_part = self.current_model.CreateModelPart("FluidMOdelPart")
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        properties = Kratos.Properties(0)
        self.fluid_model_part.AddProperties(properties)

        self.settings = Kratos.Parameters("""
        {
            "model_part_name_list" : ["FluidMOdelPart.first", "FluidMOdelPart.second"]
            "output_file_settings": {
                    "file_name"  : "test_flow_data_output.dat",
                    "folder_name": "disp_output"
            }
        }
        """)

        self.filename = "test_flow_data_output.dat"

    def testFlowsMeasuring2D_1(self):
        vel = Kratos.Array3()
        vel[0] = 0.0
        vel[1] = 1.0
        vel[2] = 0.0
        node = self.fluid_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)
        node = self.fluid_model_part.CreateNewNode(2, 1.0, 1.0, 0.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)
        node = self.fluid_model_part.CreateNewNode(3, 2.0, 2.0, 0.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)
        node = self.fluid_model_part.CreateNewNode(4, 3.0, 3.0, 0.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)

        condition_name = "LineCondition2D2N"
        self.fluid_model_part.CreateNewCondition(condition_name, 1, [1, 2], self.fluid_model_part.GetProperties()[0])
        self.fluid_model_part.CreateNewCondition(condition_name, 2, [2, 3], self.fluid_model_part.GetProperties()[0])
        self.fluid_model_part.CreateNewCondition(condition_name, 3, [3, 4], self.fluid_model_part.GetProperties()[0])

        first_smp = self.fluid_model_part.CreateSubModelPart("first")
        first_smp.AddConditions([1,2])
        second_smp = self.fluid_model_part.CreateSubModelPart("second")
        second_smp.AddConditions([3])

        flow_measurer = Fluid.FlowsMeasuringProcess(self.current_model, self.settings)
        flow_measurer.ExecuteFinalizeSolutionStep()

        with open(self.filename, "r") as f:
            f.readline()
            interesting_line = f.readline()
            splitted_line = interesting_line.split()
            self.assertEqual(splitted_line[1], "-2")
            self.assertEqual(splitted_line[2], "-1")

    def testFlowsMeasuring3D_1(self):
        vel = Kratos.Array3()
        vel[0] = 0.0
        vel[1] = 0.0
        vel[2] = 1.0
        node = self.fluid_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)
        node = self.fluid_model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)
        node = self.fluid_model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)
        node = self.fluid_model_part.CreateNewNode(4, 1.0, 1.0, 0.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)
        vel[0] = 1.0
        vel[1] = 0.0
        vel[2] = 0.0
        node = self.fluid_model_part.CreateNewNode(5, 2.0, 0.0, 0.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)
        node = self.fluid_model_part.CreateNewNode(6, 2.0, 0.0, 1.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)
        node = self.fluid_model_part.CreateNewNode(7, 2.0, 1.0, 0.0)
        node.SetSolutionStepValue(Kratos.VELOCITY, vel)

        condition_name = "SurfaceCondition3D3N"
        self.fluid_model_part.CreateNewCondition(condition_name, 1, [1, 2, 3], self.fluid_model_part.GetProperties()[0])
        self.fluid_model_part.CreateNewCondition(condition_name, 2, [3, 2, 4], self.fluid_model_part.GetProperties()[0])
        self.fluid_model_part.CreateNewCondition(condition_name, 3, [5, 6, 7], self.fluid_model_part.GetProperties()[0])

        first_smp = self.fluid_model_part.CreateSubModelPart("first")
        first_smp.AddConditions([1,2])
        second_smp = self.fluid_model_part.CreateSubModelPart("second")
        second_smp.AddConditions([3])


        flow_value_first = KratosCFD.PostProcessUtilities().ComputeFlow(first_smp)
        flow_value_second = KratosCFD.PostProcessUtilities().ComputeFlow(second_smp)

        self.assertEqual(flow_value_first, 1.0)
        self.assertEqual(flow_value_second, -0.5)

    def tearDown(self):

        try:
            os.remove(self.filename)
        except OSError:
            pass

if __name__ == '__main__':
    UnitTest.main()
