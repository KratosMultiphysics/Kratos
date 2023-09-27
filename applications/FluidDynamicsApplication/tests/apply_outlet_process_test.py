import csv

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.FluidDynamicsApplication import apply_outlet_process

class ApplyOutletProcessTest(KratosUnittest.TestCase):
    @classmethod
    def __CreateModel(self):
        model = KratosMultiphysics.Model()
        outlet_model_part = model.CreateModelPart("OutletModelPart")
        properties = outlet_model_part.CreateNewProperties(0)
        outlet_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
        outlet_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        outlet_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        outlet_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        outlet_model_part.CreateNewNode(1, 0.0, 1.0, 0.0)
        outlet_model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        outlet_model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], properties)
        for node in outlet_model_part.Nodes:
            node.AddDof(KratosMultiphysics.PRESSURE)
        return model

    @classmethod
    def __GetBlankParameters(cls):
        return KratosMultiphysics.Parameters("""{
            "name" : "Processes.KratosMultiphysics.FluidDynamicsApplication.apply_outlet_process",
            "Parameters"    : {
                "model_part_name" : "OutletModelPart"
            }
        }""")
    
    @classmethod
    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting("test_outlet_table.csv")

    def testApplyOutletProcessScalar(self):
        # Set up test
        model = self.__CreateModel()
        settings = self.__GetBlankParameters()
        
        # Customize blank settings
        settings["Parameters"].AddEmptyValue("value").SetDouble(1.0)
        
        # Create the outlet process
        outlet_process = apply_outlet_process.ApplyOutletProcess(model, settings["Parameters"])

        # Execute and check the outlet process
        self.__ExecuteAndCheckOutletProcess(model, outlet_process)

    def testApplyOutletProcessString(self):
        # Set up test
        model = self.__CreateModel()
        settings = self.__GetBlankParameters()
        
        # Customize blank settings
        settings["Parameters"].AddEmptyValue("value").SetString("1.0")
        
        # Create the outlet process
        outlet_process = apply_outlet_process.ApplyOutletProcess(model, settings["Parameters"])

        # Execute and check the outlet process
        self.__ExecuteAndCheckOutletProcess(model, outlet_process)

    def testApplyOutletProcessTable(self):
        # Set up test
        model = self.__CreateModel()
        settings = self.__GetBlankParameters()
        
        # Create an auxiliary table to be used by the process
        # This will be deleted in the tearDown function once the test is finished
        with open('test_outlet_table.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows([
                [0.0,0.0],
                [2.0,2.0]
            ])

        # Customize blank settings
        table_settings = KratosMultiphysics.Parameters("""{
            "filename" : "test_outlet_table.csv"                           
        }""")
        settings["Parameters"].AddValue("value", table_settings)
        
        # Create the outlet process
        outlet_process = apply_outlet_process.ApplyOutletProcess(model, settings["Parameters"])

        # Execute and check the outlet process
        model.GetModelPart("OutletModelPart").ProcessInfo[KratosMultiphysics.TIME] = 1.0
        self.__ExecuteAndCheckOutletProcess(model, outlet_process)

    def __ExecuteAndCheckOutletProcess(self, model, outlet_process):
        # Execute and check ExecuteInitializeSolutionStep
        outlet_process.ExecuteInitializeSolutionStep()
        for node in model.GetModelPart("OutletModelPart").Nodes:
            self.assertTrue(node.Is(KratosMultiphysics.OUTLET))
            self.assertTrue(node.IsFixed(KratosMultiphysics.PRESSURE))
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), 1.0)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE), 1.0)
        for condition in model.GetModelPart("OutletModelPart").Conditions:
            self.assertTrue(condition.Is(KratosMultiphysics.OUTLET))

        # Execute and check ExecuteFinalizeSolutionStep
        outlet_process.ExecuteFinalizeSolutionStep()
        for node in model.GetModelPart("OutletModelPart").Nodes:
            self.assertFalse(node.IsFixed(KratosMultiphysics.PRESSURE))

if __name__ == '__main__':
    KratosUnittest.main()
