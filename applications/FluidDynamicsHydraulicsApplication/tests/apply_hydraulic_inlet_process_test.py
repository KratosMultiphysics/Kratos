import csv
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities  as KratosUtilities
from KratosMultiphysics.FluidDynamicsHydraulicsApplication import apply_hydraulic_inlet_process

class ApplyHydraulicInletProcessTest(UnitTest.TestCase):

    def setUp(self):
        """Generates the model part to work on."""

        # Set up the domain model parts
        self.model = KratosMultiphysics.Model()
        test_model_part = self.model.CreateModelPart("TestModelPart")
        inlet_model_part = test_model_part.CreateSubModelPart("Inlet")
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)

        # Add VELOCITY variable to be used as DOF and DISTANCE
        test_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        test_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        
        # Modelpart: The inlet condition
        # 7 - 8 - 9
        # | \ | \ |
        # 4 - 5 - 6
        # | \ | \ |
        # 1 - 2 - 3

        # Populate model part nodes
        test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        test_model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
        test_model_part.CreateNewNode(3, 0.0, 2.0, 0.0)
        test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
        test_model_part.CreateNewNode(5, 0.0, 1.0, 1.0)
        test_model_part.CreateNewNode(6, 0.0, 2.0, 1.0)
        test_model_part.CreateNewNode(7, 0.0, 0.0, 2.0)
        test_model_part.CreateNewNode(8, 0.0, 1.0, 2.0)
        test_model_part.CreateNewNode(9, 0.0, 2.0, 2.0)
        inlet_model_part.AddNodes([1, 2, 3, 4, 5, 6, 7, 8, 9])

        # Create inlet conditions
        prop_0 = test_model_part.Properties[0]
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 4, 2], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [4, 5, 2], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 3, [2, 5, 3], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 4, [5, 6, 3], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 5, [4, 7, 5], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 6, [7, 8, 5], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 7, [6, 5, 8], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 8, [8, 9, 6], prop_0)

        # Add velocity DOFs
        for node in test_model_part.Nodes:
            node.AddDof(KratosMultiphysics.VELOCITY_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z)
            node.AddDof(KratosMultiphysics.DISTANCE)

    @classmethod
    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting("aux_table.csv")

    def testUnitHydraulicScalar3D(self):
        """Tests the process with constant discharge."""

        # Create the parabolic inlet process
        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "inlet_model_part_name" : "TestModelPart.Inlet",
                "value":1.808314132,
                "interval"        : [0.0,"End"],
                "water_depth_variable":"AUX_DISTANCE"
            }
        }""")
        process = apply_hydraulic_inlet_process.ApplyHydraulicInletProcess(self.model, settings["Parameters"])

        test_model_part = self.model.GetModelPart("TestModelPart")
        # Use inlet midpoint node as reference node
        test_node = test_model_part.GetSubModelPart("Inlet").GetNode(2)

        # Calculate normal and parallel distance
        process.ExecuteInitialize()
        self.assertEqual(test_node.Is(KratosMultiphysics.INLET), True)

        # Set current step value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(
            KratosMultiphysics.VELOCITY), [1.80831141953287, -0, -0],3)
        self.assertAlmostEqual(test_node.GetValue(KratosMultiphysics.AUX_DISTANCE), -0.5,3)

        # Remove fixity for next step
        process.ExecuteFinalizeSolutionStep()

    def testUnitHydraulicFunction3D(self):
        """Tests the process with hydrograph time dependant function."""

        # Create the parabolic inlet process
        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "inlet_model_part_name" : "TestModelPart.Inlet",
                "value":"0.5+0.2*t",
                "interval"        : [0.0,"End"],
                "water_depth_variable":"AUX_DISTANCE"
            }
        }""")
        process = apply_hydraulic_inlet_process.ApplyHydraulicInletProcess(self.model, settings["Parameters"])

        test_model_part = self.model.GetModelPart("TestModelPart")

        # Use inlet midpoint node as reference node
        test_node = test_model_part.GetSubModelPart("Inlet").GetNode(2)

        # Calculate normal and parallel distance
        process.ExecuteInitialize()
        self.assertEqual(test_node.Is(KratosMultiphysics.INLET), True)

        # Set current step value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(
            KratosMultiphysics.VELOCITY), [1.40051, 0, 0], 3)
        self.assertAlmostEqual(test_node.GetValue(
            KratosMultiphysics.AUX_DISTANCE), -0.250000875,3)

        # Remove fixity for next step
        process.ExecuteFinalizeSolutionStep()

    def testUnitHydraulicTable3D(self):
        """Tests the process with a hydrograph table."""

        # Create an auxiliary table to be used by the process
        # This will be deleted in the tearDown function once the test is finished
        with open('aux_table.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows([
                [0.0, 0.5],
                [2.0, 1.8]
            ])

        # Create the parabolic inlet process
        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "inlet_model_part_name" : "TestModelPart.Inlet",
                "value": {
                    "filename": "aux_table.csv"
                },
                "interval"        : [0.0,"End"],
                "water_depth_variable":"AUX_DISTANCE"
            }
        }""")
        process = apply_hydraulic_inlet_process.ApplyHydraulicInletProcess(self.model, settings["Parameters"])

        test_model_part = self.model.GetModelPart("TestModelPart")

        # Use inlet midpoint node as reference node
        test_node = test_model_part.GetSubModelPart("Inlet").GetNode(2)

        # Calculate normal and parallel distance
        process.ExecuteInitialize()
        self.assertEqual(test_node.Is(KratosMultiphysics.INLET), True)

        # Set current step value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(
            KratosMultiphysics.VELOCITY), [1.60761,-0,-0],3)
        self.assertAlmostEqual(test_node.GetValue(
            KratosMultiphysics.AUX_DISTANCE), -0.35742269628906254,3)

        # Remove fixity for next step
        process.ExecuteFinalizeSolutionStep()

if __name__ == '__main__':
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.PROGRESS)
    UnitTest.main()
