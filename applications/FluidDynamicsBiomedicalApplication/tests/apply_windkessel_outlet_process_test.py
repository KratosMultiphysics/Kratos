import csv

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.FluidDynamicsApplication import apply_outlet_process
from KratosMultiphysics.FluidDynamicsBiomedicalApplication import apply_windkessel_outlet_process

class ApplyWindkesselOutletProcessTest(KratosUnittest.TestCase):
    @classmethod
    def __CreateModel(self):
        model = KratosMultiphysics.Model()
        outlet_model_part = model.CreateModelPart("OutletModelPart")
        properties = outlet_model_part.CreateNewProperties(0)
        outlet_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
        outlet_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.01)
        outlet_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        outlet_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        outlet_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        outlet_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        outlet_model_part.CreateNewNode(1, 0.0, 1.0, 0.0)
        outlet_model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        outlet_model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], properties)
        for node in outlet_model_part.Nodes:
            node.AddDof(KratosMultiphysics.PRESSURE)
            node.AddDof(KratosMultiphysics.VELOCITY_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,[0.0,0.0,0.0])
        return model

    @classmethod
    def __GetBlankParameters(cls):
        return KratosMultiphysics.Parameters("""{
            "name" : "Processes.KratosMultiphysics.FluidDynamicsBiomedicalApplication.apply_windkessel_outlet_process",
            "Parameters"    : {
                "model_part_name" : "OutletModelPart"
            }
        }""")

    @classmethod
    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting("test_outlet_table.csv")

    def testApplyWindkesselOutletProcessConversionFalse(self):
        # Set up test - Check the behaviour of Windkessel update when values are written in Pa
        model = self.__CreateModel()
        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"    : "OutletModelPart",
                "mesh_id"            : 0,
                "variable_name"      : "PRESSURE",
                "constrained"        : false,
                "value"              : 0.0,
                "interval"           : [0.0,"End"],
                "resistance1"        : 0.02,
                "resistance2"        : 0.98,
                "compliance"         : 1.02,
                "venous_pressure"    : 0.0,
                "initial_pressure"   : 1.0,
                "pressure_in_mmHg"   : false,
                "echo_level"         : 1
            }
        }""")


        # Create the outlet process
        outlet_process = apply_windkessel_outlet_process.ApplyWindkesselOutletProcess(model, settings["Parameters"])

        # Execute and check the outlet process
        self.__ExecuteAndCheckOutletProcessFalse(model, outlet_process)

    def testApplyWindkesselOutletProcessConversionTrue(self):
        # Set up test - Check the behaviour of Windkessel update when values are written in mmHg and converted in Pa
        model = self.__CreateModel()
        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"    : "OutletModelPart",
                "mesh_id"            : 0,
                "variable_name"      : "PRESSURE",
                "constrained"        : false,
                "value"              : 0.0,
                "interval"           : [0.0,"End"],
                "resistance1"        : 0.02,
                "resistance2"        : 0.98,
                "compliance"         : 1.02,
                "venous_pressure"    : 0.0,
                "initial_pressure"   : 1.0,
                "pressure_in_mmHg"   : true,
                "echo_level"         : 1
            }
        }""")


        # Create the outlet process
        outlet_process = apply_windkessel_outlet_process.ApplyWindkesselOutletProcess(model, settings["Parameters"])

        # Execute and check the outlet process
        self.__ExecuteAndCheckOutletProcessTrue(model, outlet_process)

    def __ExecuteAndCheckOutletProcessFalse(self, model, outlet_process):
        # Execute and check ExecuteInitializeSolutionStep
        outlet_process.ExecuteInitializeSolutionStep()
        outlet_process.ExecuteFinalizeSolutionStep()
        for node in model.GetModelPart("OutletModelPart").Nodes:
            self.assertTrue(node.Is(KratosMultiphysics.OUTLET))
            self.assertTrue(node.IsFixed(KratosMultiphysics.PRESSURE))
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), 0.99,5)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE), 0.99,5)
        for condition in model.GetModelPart("OutletModelPart").Conditions:
            self.assertTrue(condition.Is(KratosMultiphysics.OUTLET))

    def __ExecuteAndCheckOutletProcessTrue(self, model, outlet_process):
        # Execute and check ExecuteInitializeSolutionStep
        outlet_process.ExecuteInitializeSolutionStep()
        outlet_process.ExecuteFinalizeSolutionStep()
        for node in model.GetModelPart("OutletModelPart").Nodes:
            self.assertTrue(node.Is(KratosMultiphysics.OUTLET))
            self.assertTrue(node.IsFixed(KratosMultiphysics.PRESSURE))
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), 131.5472,1)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE), 131.5472,1)
        for condition in model.GetModelPart("OutletModelPart").Conditions:
            self.assertTrue(condition.Is(KratosMultiphysics.OUTLET))

if __name__ == '__main__':
    KratosUnittest.main()
