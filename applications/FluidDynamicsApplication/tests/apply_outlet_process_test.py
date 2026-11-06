import csv

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
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
        outlet_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
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

    def testApplyOutletProcessUserDefinedCharacteristicVelocity(self):
        # Set up test with outlet velocities different from the user defined characteristic velocity
        model = self.__CreateModel()
        outlet_model_part = model.GetModelPart("OutletModelPart")
        for node in outlet_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, [float(node.Id), 0.0, 0.0])
        settings = self.__GetBlankParameters()
        settings["Parameters"].AddEmptyValue("value").SetDouble(0.0)
        settings["Parameters"].AddEmptyValue("outlet_inflow_contribution").SetBool(True)
        settings["Parameters"].AddEmptyValue("outlet_inflow_contribution_characteristic_velocity_calculation").SetString("user_defined")
        settings["Parameters"].AddEmptyValue("outlet_inflow_contribution_characteristic_velocity_value").SetDouble(2.5)

        # Check that the user defined value is kept
        outlet_process = apply_outlet_process.ApplyOutletProcess(model, settings["Parameters"])
        outlet_process.ExecuteInitializeSolutionStep()
        self.assertAlmostEqual(outlet_model_part.ProcessInfo[KratosFluid.CHARACTERISTIC_VELOCITY], 2.5)

    def testApplyOutletProcessEmptyString(self):
        # An empty function string must be rejected
        model = self.__CreateModel()
        settings = self.__GetBlankParameters()
        settings["Parameters"].AddEmptyValue("value").SetString("")
        with self.assertRaisesRegex(Exception, "function string is empty"):
            apply_outlet_process.ApplyOutletProcess(model, settings["Parameters"])

    def testApplyOutletProcessHydrostaticPropertiesDensity(self):
        # Hydrostatic outlet without DENSITY in the nodal data (e.g. QSVMS): the density is taken from the properties
        model = self.__CreateHydrostaticModel(nodal_density = False, densities = [1000.0])
        outlet_process = apply_outlet_process.ApplyOutletProcess(model, self.__GetHydrostaticParameters()["Parameters"])
        outlet_process.ExecuteInitializeSolutionStep()
        self.__CheckHydrostaticPressure(model, 1000.0)

    def testApplyOutletProcessHydrostaticNodalDensity(self):
        # Hydrostatic outlet with DENSITY in the nodal data: the nodal value is used
        model = self.__CreateHydrostaticModel(nodal_density = True, densities = [1000.0])
        outlet_process = apply_outlet_process.ApplyOutletProcess(model, self.__GetHydrostaticParameters()["Parameters"])
        outlet_process.ExecuteInitializeSolutionStep()
        self.__CheckHydrostaticPressure(model, 1.0)

    def testApplyOutletProcessHydrostaticNonUniqueDensity(self):
        # Hydrostatic outlet without DENSITY in the nodal data and different densities in the properties
        model = self.__CreateHydrostaticModel(nodal_density = False, densities = [1000.0, 1.2])
        outlet_process = apply_outlet_process.ApplyOutletProcess(model, self.__GetHydrostaticParameters()["Parameters"])
        with self.assertRaisesRegex(Exception, "unique DENSITY"):
            outlet_process.ExecuteInitializeSolutionStep()

    @classmethod
    def __CreateHydrostaticModel(cls, nodal_density, densities):
        # Outlet with the nodes at y = 1 (top) and y = 0 (bottom) and the body force along -y
        model = KratosMultiphysics.Model()
        outlet_model_part = model.CreateModelPart("OutletModelPart")
        outlet_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
        for variable in [KratosMultiphysics.PRESSURE, KratosMultiphysics.EXTERNAL_PRESSURE, KratosMultiphysics.BODY_FORCE]:
            outlet_model_part.AddNodalSolutionStepVariable(variable)
        if nodal_density:
            outlet_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        for i, density in enumerate(densities):
            outlet_model_part.CreateNewProperties(i + 1).SetValue(KratosMultiphysics.DENSITY, density)
        outlet_model_part.CreateNewNode(1, 0.0, 1.0, 0.0)
        outlet_model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        outlet_model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], outlet_model_part.CreateNewProperties(0))
        for node in outlet_model_part.Nodes:
            node.AddDof(KratosMultiphysics.PRESSURE)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE, [0.0, -9.81, 0.0])
            if nodal_density:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 1.0)
        return model

    @classmethod
    def __GetHydrostaticParameters(cls):
        return KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name" : "OutletModelPart",
                "value" : 0.0,
                "hydrostatic_outlet" : true,
                "h_top" : 0.0
            }
        }""")

    def __CheckHydrostaticPressure(self, model, density):
        # Hydrostatic pressure: zero at the top node and density * g * 1.0 at the bottom one (relative tolerance)
        expected = {1: 0.0, 2: density * 9.81 * 1.0}
        for node in model.GetModelPart("OutletModelPart").Nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE), expected[node.Id], delta = 1.0e-9 * max(1.0, expected[node.Id]))
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE), expected[node.Id], delta = 1.0e-9 * max(1.0, expected[node.Id]))

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
