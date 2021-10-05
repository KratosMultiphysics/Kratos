import KratosMultiphysics
from KratosMultiphysics import FluidDynamicsApplication
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.FluidDynamicsApplication import apply_compressible_navier_stokes_boundary_conditions_process


class ApplyMachDependentBoundaryConditionsTest(UnitTest.TestCase):

    def setUp(self):
        """Generates the model part to work on."""
        self.model = KratosMultiphysics.Model()

        main_mpart = self.model.CreateModelPart("main_model_part")
        open_mpart = main_mpart.CreateSubModelPart("open_boundaries")
        closed_mpart = main_mpart.CreateSubModelPart("closed_boundaries")
        inlet_mpart = open_mpart.CreateSubModelPart("inlet")
        outlet_mpart = open_mpart.CreateSubModelPart("outlet")

        main_mpart.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        main_mpart.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        main_mpart.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        main_mpart.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        main_mpart.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM)

        # Modelpart:
        # 7 - 8 - 9  : outlet
        # |   |   |       Δ
        # 4 - 5 - 6       |  direction of flow
        # |   |   |       |
        # 1 - 2 - 3  :  inlet

        main_mpart.CreateNewNode(1, 0.0, 0.0, 0.0)
        main_mpart.CreateNewNode(2, 1.0, 0.0, 0.0)
        main_mpart.CreateNewNode(3, 2.0, 0.0, 0.0)
        main_mpart.CreateNewNode(4, 0.0, 1.0, 0.0)
        main_mpart.CreateNewNode(5, 1.0, 1.0, 0.0)
        main_mpart.CreateNewNode(6, 2.0, 1.0, 0.0)
        main_mpart.CreateNewNode(7, 0.0, 2.0, 0.0)
        main_mpart.CreateNewNode(8, 1.0, 2.0, 0.0)
        main_mpart.CreateNewNode(9, 2.0, 2.0, 0.0)

        inlet_mpart.AddNodes([1, 2, 3])
        outlet_mpart.AddNodes([7, 8, 9])
        closed_mpart.AddNodes([4, 6])

        props = main_mpart.Properties[0]

        inlet_mpart.CreateNewCondition("LineCondition2D2N", 1, [1, 2], props)
        inlet_mpart.CreateNewCondition("LineCondition2D2N", 2, [2, 3], props)

        main_mpart.CreateNewCondition("LineCondition2D2N", 3, [3, 6], props)
        main_mpart.CreateNewCondition("LineCondition2D2N", 4, [6, 9], props)

        outlet_mpart.CreateNewCondition("LineCondition2D2N", 5, [9, 8], props)
        outlet_mpart.CreateNewCondition("LineCondition2D2N", 6, [8, 7], props)

        main_mpart.CreateNewCondition("LineCondition2D2N", 7, [7, 4], props)
        main_mpart.CreateNewCondition("LineCondition2D2N", 8, [4, 1], props)

        for node in main_mpart.Nodes:
            node.AddDof(KratosMultiphysics.TEMPERATURE)
            node.AddDof(KratosMultiphysics.DENSITY)

            node.AddDof(KratosMultiphysics.VELOCITY_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z)

        self._ResetVariablesToZero()

    def _ResetVariablesToZero(self):
        for node in self.model["main_model_part"].Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0.0)

            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, 0.0)

            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_Y, 1.0)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_Z, 0.0)

    def testDoubleVariable(self):
        """Tests that the process works well with variables of type array_1d<double, 3> works well with variables of type double."""
        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "model_part_name" : "main_model_part",
                "flow_direction_variable" : "MOMENTUM",
                "subsonic_boundary_conditions" : [
                    {
                        "variable_name" : "DENSITY",
                        "value" : 1.225,
                        "interval" : [0, "End"]
                    }
                ],
                "supersonic_boundary_conditions" : [
                    {
                        "variable_name" : "TEMPERATURE",
                        "value" : 273.15,
                        "interval" : [0, "End"]
                    }
                ]
            }
        }
        """)

        process = apply_compressible_navier_stokes_boundary_conditions_process.Factory(settings, self.model)
        main_model_part = self.model["main_model_part"]

        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()

        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)

        # Subsonic case
        for node in main_model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 0.0)

        process.ExecuteInitializeSolutionStep()

        for node in main_model_part.Nodes:
            self.assertTrue(node.IsFixed(KratosMultiphysics.DENSITY),
                msg="Failed to fix subsonic boundary condition (Node #{}).".format(node.Id))
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), 1.225,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

            self.assertFalse(node.IsFixed(KratosMultiphysics.TEMPERATURE),
                msg="Mistakenly fixed supersonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 0.0,
                msg="Mistakenly modified value for supersonic boundary condition (Node #%d)." % node.Id)

        self._ResetVariablesToZero()
        process.ExecuteFinalizeSolutionStep()

        # Supersonic case
        for node in main_model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 5.0)

        process.ExecuteInitializeSolutionStep()

        for node in self.model["main_model_part.open_boundaries"].Nodes:
            self.assertTrue(node.IsFixed(KratosMultiphysics.TEMPERATURE),
                msg="Failed to fix supersonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 273.15,
                msg="Failed to set value for supersonic boundary condition (Node #%d)." % node.Id)

            self.assertFalse(node.IsFixed(KratosMultiphysics.DENSITY),
                msg="Mistakenly fixed subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), 0.0,
                msg="Mistakenly modified value for subsonic boundary condition (Node #%d)." % node.Id)

        for node in self.model["main_model_part.closed_boundaries"].Nodes:
            # These boundaries are still subsonic because they're perpendicular to the flow
            self.assertTrue(node.IsFixed(KratosMultiphysics.DENSITY),
                msg="Failed to fix subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), 1.225,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

            self.assertFalse(node.IsFixed(KratosMultiphysics.TEMPERATURE),
                msg="Mistakenly fixed supersonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 0.0,
                msg="Mistakenly modified value for supersonic boundary condition (Node #%d)." % node.Id)

        process.ExecuteFinalizeSolutionStep()

    def testVectorVariable(self):
        """Tests that the process works well with variables of type array_1d<double, 3>."""
        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "model_part_name" : "main_model_part",
                "flow_direction_variable" : "MOMENTUM",
                "subsonic_boundary_conditions" : [
                    {
                        "variable_name" : "VELOCITY",
                        "value" : [1.0, 2.0, 3.0],
                        "constrained" : [false, false, true],
                        "interval" : [0, "End"]
                    }
                ],
                "supersonic_boundary_conditions" : [
                    {
                        "variable_name" : "TEMPERATURE",
                        "value" : 273.15,
                        "interval" : [0, "End"]
                    }
                ]
            }
        }
        """)

        process = apply_compressible_navier_stokes_boundary_conditions_process.Factory(settings, self.model)
        model_part = self.model["main_model_part"]

        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()

        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)

        # Subsonic case
        for node in model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 0.5)

        process.ExecuteInitializeSolutionStep()

        for node in self.model["main_model_part.open_boundaries"].Nodes:
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_X),
                msg="Mistakenly fixed unconstrained component of vector variable (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable (Node #%d)." % node.Id)

            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_Y),
                msg="Mistakenly fixed unconstrained component of vector variable (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable (Node #%d)." % node.Id)

            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_Z),
                msg="Failed to fix subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z), 3.0,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

            self.assertFalse(node.IsFixed(KratosMultiphysics.TEMPERATURE),
                msg="Mistakenly fixed supersonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 0.0,
                msg="Mistakenly modified value for supersonic boundary condition (Node #%d)." % node.Id)

        self._ResetVariablesToZero()
        process.ExecuteFinalizeSolutionStep()

        # Supersonic case
        for node in model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 15)

        process.ExecuteInitializeSolutionStep()

        for node in self.model["main_model_part.open_boundaries"].Nodes:
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_X),
                msg="Mistakenly fixed unconstrained component of vector variable (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable (Node #%d)." % node.Id)

            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_Y),
                msg="Mistakenly fixed unconstrained component of vector variable (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable (Node #%d)." % node.Id)

            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_Z),
                msg="Mistakenly fixed subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z), 0.0,
                msg="Mistakenly modified value for subsonic boundary condition (Node #%d)." % node.Id)

            self.assertTrue(node.IsFixed(KratosMultiphysics.TEMPERATURE),
                msg="Failed to fix supersonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 273.15,
                msg="Failed to set value for supersonic boundary condition (Node #%d)." % node.Id)

        for node in self.model["main_model_part.closed_boundaries"].Nodes:
            # These boundaries are still subsonic because they're perpendicular to the flow
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_X),
                msg="Mistakenly fixed unconstrained component of vector variable (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable (Node #%d)." % node.Id)

            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_Y),
                msg="Mistakenly fixed unconstrained component of vector variable (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable (Node #%d)." % node.Id)

            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_Z),
                msg="Failed to fix subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z), 3.0,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

            self.assertFalse(node.IsFixed(KratosMultiphysics.TEMPERATURE),
                msg="Mistakenly fixed supersonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 0.0,
                msg="Mistakenly modified value for supersonic boundary condition (Node #%d)." % node.Id)

        process.ExecuteFinalizeSolutionStep()

    def testInterval(self):
        """Tests that the boundary conditions are enforced only during the relevant time interval."""
        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "model_part_name" : "main_model_part",
                "flow_direction_variable" : "MOMENTUM",
                "subsonic_boundary_conditions" : [
                    {
                        "variable_name" : "DENSITY",
                        "value" : 1.225,
                        "interval" : [0, 5.0]
                    }
                ]
            }
        }
        """)

        process = apply_compressible_navier_stokes_boundary_conditions_process.Factory(settings, self.model)
        model_part = self.model["main_model_part"]

        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()

        for node in model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 0.5)

        # t inside interval
        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)


        process.ExecuteInitializeSolutionStep()

        for node in model_part.Nodes:
            self.assertTrue(node.IsFixed(KratosMultiphysics.DENSITY),
                msg="Failed to fix subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), 1.225,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

        self._ResetVariablesToZero()
        process.ExecuteFinalizeSolutionStep()

        # t outside interval
        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 100.0)

        process.ExecuteInitializeSolutionStep()

        for node in model_part.Nodes:
            self.assertFalse(node.IsFixed(KratosMultiphysics.DENSITY),
                msg="Mistakenly fixed boundary condition outside the specified time interval (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), 0.0,
                msg="Mistakenly modified value for boundary condition outside the specified time interval (Node #%d)." % node.Id)

        process.ExecuteFinalizeSolutionStep()

    def testErrorMissingValue(self):
        """Ensures an easy-to understand error is raised when the parameter 'value' is missing."""
        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "model_part_name" : "main_model_part",
                "flow_direction_variable" : "MOMENTUM",
                "subsonic_boundary_conditions" : [
                    {
                        "variable_name" : "DENSITY",
                        "interval" : [0, 5.0]
                    }
                ]
            }
        }
        """)

        with self.assertRaises(RuntimeError) as context:
            _ = apply_compressible_navier_stokes_boundary_conditions_process.Factory(settings, self.model)

        self.assertIn("Missing double (or vector) field 'value' in boundary condition", str(context.exception))
        self.assertIn("DENSITY", str(context.exception))

    def testErrorMissingConstraint(self):
        """Ensures an easy-to understand error is raised when the parameter 'constraint' is missing."""
        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "model_part_name" : "main_model_part",
                "flow_direction_variable" : "MOMENTUM",
                "subsonic_boundary_conditions" : [
                    {
                        "variable_name" : "VELOCITY",
                        "value" : [0,1,3],
                        "interval" : [0, 5.0]
                    }
                ]
            }
        }
        """)

        with self.assertRaises(RuntimeError) as context:
            _ = apply_compressible_navier_stokes_boundary_conditions_process.Factory(settings, self.model)

        self.assertIn("Missing boolean array field \'constrained", str(context.exception))
        self.assertIn("VELOCITY", str(context.exception))

    def testErrorTooManyValues(self):
        """Raising this exception is particularly important because it prevents possible segmentation faults."""
        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "model_part_name" : "main_model_part",
                "flow_direction_variable" : "MOMENTUM",
                "subsonic_boundary_conditions" : [
                    {
                        "variable_name" : "VELOCITY",
                        "value" : [0,1, 3, 4],
                        "constrained" : [true, true, true, true],
                        "interval" : [0, 5.0]
                    }
                ]
            }
        }
        """)

        with self.assertRaises(RuntimeError) as context:
            _ = apply_compressible_navier_stokes_boundary_conditions_process.Factory(settings, self.model)

        self.assertIn("Allowed vector variables are at most 3-dimensional", str(context.exception))
        self.assertIn("VELOCITY", str(context.exception))


    def testSameVariable(self):
        """Ensures fixing diferent axes on the same variable does not cause unexpected interfearances."""
        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "model_part_name" : "main_model_part",
                "flow_direction_variable" : "MOMENTUM",
                "subsonic_boundary_conditions" : [
                    {
                        "variable_name" : "VELOCITY",
                        "value" : [1.0, 2.0, 3.0],
                        "constrained" : [true, true, true],
                        "interval" : [0, "End"]
                    }
                ],
                "supersonic_boundary_conditions" : [
                    {
                        "variable_name" : "VELOCITY",
                        "value" : [-1.0, -2.0, -3.0],
                        "constrained" : [false, false, true],
                        "interval" : [0, "End"]
                    }
                ]
            }
        }
        """)

        process = apply_compressible_navier_stokes_boundary_conditions_process.Factory(settings, self.model)
        main_model_part = self.model["main_model_part"]

        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()

        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)

        # Subsonic case
        for node in main_model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 0.0)

        process.ExecuteInitializeSolutionStep()

        for node in self.model["main_model_part.open_boundaries"].Nodes:
            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_X),
                msg="Failed to fix subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X), 1.0,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_Y),
                msg="Failed to fix subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y), 2.0,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_Z),
                msg="Failed to fix subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z), 3.0,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

        self._ResetVariablesToZero()
        process.ExecuteFinalizeSolutionStep()

        # Supersonic case
        for node in main_model_part.Nodes:
            node.SetValue(FluidDynamicsApplication.MACH, 15)

        process.ExecuteInitializeSolutionStep()

        for node in self.model["main_model_part.open_boundaries"].Nodes:
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_X),
                msg="Mistakenly fixed unconstrained component of vector variable (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable (Node #%d)." % node.Id)

            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_Y),
                msg="Mistakenly fixed unconstrained component of vector variable (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y), 0.0,
                msg="Mistakenly modified value for unconstrained component of vector variable (Node #%d)." % node.Id)

            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_Z),
                msg="Failed to fix supersonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z), -3.0,
                msg="Failed to set value for supersonic boundary condition (Node #%d)." % node.Id)

        for node in self.model["main_model_part.closed_boundaries"].Nodes:
            # These boundaries are still subsonic because they're perpendicular to the flow
            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_X),
                msg="Failed to fix subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X), 1.0,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_Y),
                msg="Failed to fix subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y), 2.0,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_Z),
                msg="Failed to fix subsonic boundary condition (Node #%d)." % node.Id)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z), 3.0,
                msg="Failed to set value for subsonic boundary condition (Node #%d)." % node.Id)

        process.ExecuteFinalizeSolutionStep()


if __name__ == '__main__':
    UnitTest.main()
