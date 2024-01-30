import csv

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities  as KratosUtilities
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.FluidDynamicsBiomedicalApplication as KratosBio
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid
from KratosMultiphysics.FluidDynamicsBiomedicalApplication import apply_parabolic_inlet_process

class ApplyParabolicInletProcessTest(UnitTest.TestCase):

    def setUp(self):
        """Generates the model part to work on."""

        # Set up the domain model parts
        self.model = KratosMultiphysics.Model()
        test_model_part = self.model.CreateModelPart("TestModelPart")
        inlet_model_part = test_model_part.CreateSubModelPart("Inlet")
        outlet_model_part = test_model_part.CreateSubModelPart("Outlet")
        wall_model_part = test_model_part.CreateSubModelPart("Wall")

        # Add VELOCITY variable to be used as DOF
        test_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        # Modelpart:
        # 7 - 8 - 9  : outlet
        # |   |   |       Δ
        # 4 - 5 - 6       |  direction of flow
        # |   |   |       |
        # 1 - 2 - 3  :  inlet

        # Populate model part nodes
        test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        test_model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
        test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        test_model_part.CreateNewNode(5, 1.0, 1.0, 0.0)
        test_model_part.CreateNewNode(6, 2.0, 1.0, 0.0)
        test_model_part.CreateNewNode(7, 0.0, 2.0, 0.0)
        test_model_part.CreateNewNode(8, 1.0, 2.0, 0.0)
        test_model_part.CreateNewNode(9, 2.0, 2.0, 0.0)

        inlet_model_part.AddNodes([1, 2, 3])
        outlet_model_part.AddNodes([7, 8, 9])
        wall_model_part.AddNodes([1, 3, 4, 6, 7, 9])

        # Create wall conditions
        prop_0 = test_model_part.Properties[0]
        wall_model_part.CreateNewCondition("LineCondition2D2N", 1, [4, 1], prop_0)
        wall_model_part.CreateNewCondition("LineCondition2D2N", 3, [3, 6], prop_0)
        wall_model_part.CreateNewCondition("LineCondition2D2N", 2, [7, 4], prop_0)
        wall_model_part.CreateNewCondition("LineCondition2D2N", 4, [6, 9], prop_0)

        # Create inlet conditions
        inlet_model_part.CreateNewCondition("LineCondition2D2N", 5, [1, 2], prop_0)
        inlet_model_part.CreateNewCondition("LineCondition2D2N", 6, [2, 3], prop_0)

        # Create outlet conditions
        outlet_model_part.CreateNewCondition("LineCondition2D2N", 7, [8, 7], prop_0)
        outlet_model_part.CreateNewCondition("LineCondition2D2N", 8, [9, 8], prop_0)

        # Create elements
        prop_1 = test_model_part.Properties[1]
        test_model_part.CreateNewElement("Element2D3N", 1, [1,2,4], prop_1)
        test_model_part.CreateNewElement("Element2D3N", 2, [2,5,4], prop_1)
        test_model_part.CreateNewElement("Element2D3N", 3, [2,3,5], prop_1)
        test_model_part.CreateNewElement("Element2D3N", 4, [3,6,5], prop_1)
        test_model_part.CreateNewElement("Element2D3N", 5, [4,5,7], prop_1)
        test_model_part.CreateNewElement("Element2D3N", 6, [5,8,7], prop_1)
        test_model_part.CreateNewElement("Element2D3N", 7, [5,6,8], prop_1)
        test_model_part.CreateNewElement("Element2D3N", 8, [6,9,8], prop_1)

        # Add velocity DOFs
        for node in test_model_part.Nodes:
            node.AddDof(KratosMultiphysics.VELOCITY_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z)

        # Call the fluid check and prepare model to calculate the neigbours
        settings = KratosMultiphysics.Parameters("""{
            "volume_model_part_name" : "TestModelPart",
            "skin_parts" : ["Wall", "Inlet", "Outlet"],
            "assign_neighbour_elements_to_conditions" : true
        }""")
        check_and_prepare_proc = check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(
            test_model_part,
            settings)
        check_and_prepare_proc.Execute()

    @classmethod
    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting("aux_table.csv")

    def testUnitParabolaScalar2D(self):
        """Tests the process with a unit maximum value parabola."""

        # Create the parabolic inlet process
        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "wall_model_part_name": "TestModelPart.Wall",
                "inlet_model_part_name": "TestModelPart.Inlet",
                "value" : 1.0,
                "value_is_average" : false,
                "value_is_flow_rate" : false
            }
        }""")
        process = apply_parabolic_inlet_process.Factory(settings, self.model)

        test_model_part = self.model.GetModelPart("TestModelPart")
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        # Use inlet midpoint node as reference node
        test_node = test_model_part.GetSubModelPart("Inlet").GetNode(2)

        # Calculate normal and parallel distance
        process.ExecuteInitialize()
        self.assertEqual(test_node.Is(KratosMultiphysics.INLET), True)
        self.assertVectorAlmostEqual(test_node.GetValue(KratosCFD.INLET_NORMAL), [0.0,-1.0,0.0])
        self.assertLessEqual(abs(test_node.GetValue(KratosBio.WALL_DISTANCE) - 0.9999971715808751), 1.0e-6)

        # Set initial value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        process.ExecuteBeforeSolutionLoop()
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.0,1.0,0.0])

        # Set current step value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.0,1.0,0.0])

        # Remove fixity for next step
        process.ExecuteFinalizeSolutionStep()
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))

    def testUnitParabolaAverage2D(self):
        """Tests the process with a unit maximum value parabola."""

        # Create the parabolic inlet process
        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "wall_model_part_name": "TestModelPart.Wall",
                "inlet_model_part_name": "TestModelPart.Inlet",
                "value" : 0.5,
                "value_is_average" : true,
                "value_is_flow_rate" : false
            }
        }""")
        process = apply_parabolic_inlet_process.Factory(settings, self.model)

        test_model_part = self.model.GetModelPart("TestModelPart")
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        # Use inlet midpoint node as reference node
        test_node = test_model_part.GetSubModelPart("Inlet").GetNode(2)

        # Calculate normal and parallel distance
        process.ExecuteInitialize()
        self.assertEqual(test_node.Is(KratosMultiphysics.INLET), True)
        self.assertVectorAlmostEqual(test_node.GetValue(KratosCFD.INLET_NORMAL), [0.0,-1.0,0.0])
        self.assertLessEqual(abs(test_node.GetValue(KratosBio.WALL_DISTANCE) - 0.9999971715808751), 1.0e-6)

        # Set initial value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        process.ExecuteBeforeSolutionLoop()
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.0,1.0,0.0])

        # Set current step value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.0,1.0,0.0])

        # Remove fixity for next step
        process.ExecuteFinalizeSolutionStep()
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))

    def testUnitParabolaString2D(self):
        """Tests the process with a string maximum value parabola."""

        # Create the parabolic inlet process
        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "wall_model_part_name": "TestModelPart.Wall",
                "inlet_model_part_name": "TestModelPart.Inlet",
                "value" : "2.0*t*x",
                "value_is_average" : false,
                "value_is_flow_rate" : false
            }
        }""")
        process = apply_parabolic_inlet_process.Factory(settings, self.model)

        test_model_part = self.model.GetModelPart("TestModelPart")
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        # Use inlet midpoint node as reference node
        test_node = test_model_part.GetSubModelPart("Inlet").GetNode(2)

        # Calculate normal and parallel distance
        process.ExecuteInitialize()
        self.assertEqual(test_node.Is(KratosMultiphysics.INLET), True)
        self.assertVectorAlmostEqual(test_node.GetValue(KratosCFD.INLET_NORMAL), [0.0,-1.0,0.0])
        self.assertLessEqual(abs(test_node.GetValue(KratosBio.WALL_DISTANCE) - 0.9999971715808751), 1.0e-6)

        # Set initial value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        process.ExecuteBeforeSolutionLoop()
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.0,0.0,0.0])

        # Set current step value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.0,2.0,0.0])

        # Remove fixity for next step
        process.ExecuteFinalizeSolutionStep()
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))

    def testUnitParabolaFlowRateTable2D(self):
        """Tests the process with a table maximum flow rate value parabola."""

        # Create an auxiliary table to be used by the process
        # This will be deleted in the tearDown function once the test is finished
        with open('aux_table.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows([
                [0.0,0.0],
                [2.0,1.0]
            ])

        # Create the parabolic inlet process
        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "wall_model_part_name": "TestModelPart.Wall",
                "inlet_model_part_name": "TestModelPart.Inlet",
                "value" : {
                    "filename" : "aux_table.csv"
                },
                "value_is_average" : false,
                "value_is_flow_rate" : true
            }
        }""")
        process = apply_parabolic_inlet_process.Factory(settings, self.model)

        test_model_part = self.model.GetModelPart("TestModelPart")
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        # Use inlet midpoint node as reference node
        test_node = test_model_part.GetSubModelPart("Inlet").GetNode(2)

        # Calculate normal and parallel distance
        process.ExecuteInitialize()
        self.assertEqual(test_node.Is(KratosMultiphysics.INLET), True)
        self.assertVectorAlmostEqual(test_node.GetValue(KratosCFD.INLET_NORMAL), [0.0,-1.0,0.0])
        self.assertLessEqual(abs(test_node.GetValue(KratosBio.WALL_DISTANCE) - 0.9999971715808751), 1.0e-6)

        # Set initial value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        process.ExecuteBeforeSolutionLoop()
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.0,0.0,0.0])

        # Set current step value
        test_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertTrue(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))
        self.assertVectorAlmostEqual(test_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.0,0.25,0.0])

        # Remove fixity for next step
        process.ExecuteFinalizeSolutionStep()
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_X))
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_Y))
        self.assertFalse(test_node.IsFixed(KratosMultiphysics.VELOCITY_Z))

if __name__ == '__main__':
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.PROGRESS)
    UnitTest.main()
