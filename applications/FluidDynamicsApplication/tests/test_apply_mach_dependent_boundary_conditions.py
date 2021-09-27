import KratosMultiphysics
from KratosMultiphysics import FluidDynamicsApplication
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.FluidDynamicsApplication import apply_mach_depenedent_boundary_conditions

class ApplyMachDependentBoundaryConditionsTest(UnitTest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        mpart = self.model.CreateModelPart("main_model_part")

        mpart.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        mpart.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        
        mpart.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        mpart.CreateNewNode(1, 0.0, 0.0, 0.0)
        mpart.CreateNewNode(2, 1.0, 0.0, 0.0)   
        mpart.CreateNewNode(3, 0.0, 1.0, 0.0)
        mpart.CreateNewNode(4, 0.0, 0.0, 1.0)

        for node in mpart.Nodes:
            node.AddDof(KratosMultiphysics.TEMPERATURE)
            node.AddDof(KratosMultiphysics.DENSITY)

            node.AddDof(KratosMultiphysics.VELOCITY_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z)

            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, 0.0)
        
        self._ResetVariablesToZero()

    def _ResetVariablesToZero(self):
        for node in self.model["main_model_part"].Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, 0.0)

    def testDoubleVariable(self):
        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "model_part_name" : "main_model_part",
                "subsonic_boundary_conditions" : [
                    {
                        "variable" : "DENSITY",
                        "value" : 1.225,
                        "interval" : [0, "End"]
                    }
                ],
                "supersonic_boundary_conditions" : [
                    {
                        "variable" : "TEMPERATURE",
                        "value" : 273.15,
                        "interval" : [0, "End"]
                    }
                ]
            }
        }
        """)

        process = apply_mach_depenedent_boundary_conditions.Factory(settings, self.model)
        model_part = self.model["main_model_part"]

        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()

        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)

        # Subsonic case
        for node in model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 0.5)

        process.ExecuteInitializeSolutionStep()

        for node in model_part.Nodes:
            self.assertTrue(node.IsFixed(KratosMultiphysics.DENSITY),
                msg="Failed to fix subsonic boundary condition.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), 1.225,
                msg="Failed to set value for subsonic boundary condition.")

            self.assertFalse(node.IsFixed(KratosMultiphysics.TEMPERATURE),
                msg="Mistakenly fixed supersonic boundary condition.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 0.0,
                msg="Mistakenly modified value for supersonic boundary condition.")

        self._ResetVariablesToZero()
        process.ExecuteFinalizeSolutionStep()

        # Supersonic case
        for node in model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 1.5)

        process.ExecuteInitializeSolutionStep()

        for node in model_part.Nodes:
            self.assertTrue(node.IsFixed(KratosMultiphysics.TEMPERATURE),
                msg="Failed to fix supersonic boundary condition.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 273.15,
                msg="Failed to set value for supersonic boundary condition.")

            self.assertFalse(node.IsFixed(KratosMultiphysics.DENSITY),
                msg="Mistakenly fixed subsonic boundary condition.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), 0.0,
                msg="Mistakenly modified value for subsonic boundary condition.")
        
        process.ExecuteFinalizeSolutionStep()

    def testVectorVariable(self):
        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "model_part_name" : "main_model_part",
                "subsonic_boundary_conditions" : [
                    {
                        "variable" : "VELOCITY",
                        "value" : [1.0, 2.0, 3.0],
                        "constrained" : [false, false, true],
                        "interval" : [0, "End"]
                    }
                ],
                "supersonic_boundary_conditions" : [
                    {
                        "variable" : "TEMPERATURE",
                        "value" : 273.15,
                        "interval" : [0, "End"]
                    }
                ]
            }
        }
        """)

        process = apply_mach_depenedent_boundary_conditions.Factory(settings, self.model)
        model_part = self.model["main_model_part"]

        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()

        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)

        # Subsonic case
        for node in model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 0.5)

        process.ExecuteInitializeSolutionStep()

        for node in model_part.Nodes:
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_X),
                msg="Mistakenly fixed unconstrained component of vector variable.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable.")

            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_Y),
                msg="Mistakenly fixed unconstrained component of vector variable.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable.")

            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_Z),
                msg="Failed to fix subsonic boundary condition.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z), 3.0,
                msg="Failed to set value for subsonic boundary condition.")

            self.assertFalse(node.IsFixed(KratosMultiphysics.TEMPERATURE),
                msg="Mistakenly fixed supersonic boundary condition.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 0.0,
                msg="Mistakenly modified value for supersonic boundary condition.")

        self._ResetVariablesToZero()
        process.ExecuteFinalizeSolutionStep()

        # Supersonic case
        for node in model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 1.5)

        process.ExecuteInitializeSolutionStep()

        for node in model_part.Nodes:
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_X),
                msg="Mistakenly fixed unconstrained component of vector variable.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable.")

            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_Y),
                msg="Mistakenly fixed unconstrained component of vector variable.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable.")

            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_Z),
                msg="Mistakenly fixed subsonic boundary condition.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z), 0.0,
                msg="Mistakenly modified value for subsonic boundary condition.")

            self.assertTrue(node.IsFixed(KratosMultiphysics.TEMPERATURE),
                msg="Failed to fix supersonic boundary condition.")
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 273.15,
                msg="Failed to set value for supersonic boundary condition.")
        
        process.ExecuteFinalizeSolutionStep()


if __name__ == '__main__':
    UnitTest.main()
