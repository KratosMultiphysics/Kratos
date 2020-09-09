import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as Fluid
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as k_utils
import KratosMultiphysics.FluidDynamicsApplication.flow_output_process as flow_process

class FlowsMeasuringUtilityTest(UnitTest.TestCase):

    def setUp(self):
        self.current_model = Kratos.Model()
        self.fluid_model_part = self.current_model.CreateModelPart("FluidMOdelPart")
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        properties = Kratos.Properties(0)
        self.fluid_model_part.AddProperties(properties)

        self.settings = Kratos.Parameters("""
        {
            "model_part_name_list" : ["FluidMOdelPart.first", "FluidMOdelPart.second"],
            "output_file_settings": {
                    "file_name"  : "test_flow_data_output.dat",
                    "folder_name": ""
            }
        }
        """)

        self.filename = self.settings["output_file_settings"]["file_name"].GetString()
        print(self.filename)

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

        flow_measurer = flow_process.FlowOutputProcess(self.current_model, self.settings)
        flow_measurer.ExecuteInitialize()
        flow_measurer.ExecuteInitializeSolutionStep()
        flow_measurer.ExecuteFinalizeSolutionStep()
        flow_measurer.ExecuteFinalize()

        with open(self.filename, "r") as f:
            lines=f.readlines()
            interesting_line = lines[2]
            splitted_line = interesting_line.split()
            self.assertAlmostEqual(float(splitted_line[1]), -2.0)
            self.assertAlmostEqual(float(splitted_line[2]), -1.0)

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


        flow_value_first = Fluid.FluidPostProcessUtilities().CalculateFlow(first_smp)
        flow_value_second = Fluid.FluidPostProcessUtilities().CalculateFlow(second_smp)

        self.assertAlmostEqual(flow_value_first, 1.0)
        self.assertAlmostEqual(flow_value_second, -0.5)

    def tearDown(self):
        k_utils.DeleteFileIfExisting(self.filename)


if __name__ == '__main__':
    UnitTest.main()
