import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication import apply_mass_flow_process

class ApplyMassFlowProcessTest(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("FluidModelPart")
        self.model_part.SetBufferSize(2)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        with KratosUnittest.WorkFolderScope(".", __file__):
            KratosMultiphysics.ModelPartIO("Cavity/square5").ReadModelPart(self.model_part)
        for var in [KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.VELOCITY_Z]:
            KratosMultiphysics.VariableUtils().AddDof(var, self.model_part)

        # Orient the conditions outwards as the fluid solvers do
        tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        tmoc(self.model_part, throw_errors, flags).Execute()

        # Unit square with the inlet in the left wall and the outlet in the right one
        self.inlet_model_part = self.model_part.GetSubModelPart("NoSlip2D_left_wall")
        self.outlet_model_part = self.model_part.GetSubModelPart("NoSlip2D_right_wall")

        self.density = 1.5

    def test_uniform(self):
        process = self.__CreateProcess()
        process.ExecuteInitializeSolutionStep()

        self.assertAlmostEqual(self.__CalculateMassFlow(self.inlet_model_part), 2.0, 12)
        self.assertAlmostEqual(self.__CalculateMassFlow(self.outlet_model_part), -2.0, 12)
        for node in list(self.inlet_model_part.Nodes) + list(self.outlet_model_part.Nodes):
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [2.0 / self.density, 0.0, 0.0], 12)
            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_X))
            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_Y))

        process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_X))
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_Y))

    def test_uniform_not_constrained(self):
        process = self.__CreateProcess(constrained=False)
        process.ExecuteInitializeSolutionStep()

        self.assertAlmostEqual(self.__CalculateMassFlow(self.inlet_model_part), 2.0, 12)
        for node in self.model_part.Nodes:
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_X))

    def test_uniform_respects_fixed_nodes(self):
        # Emulate the no-slip processes of the top and bottom walls (these share the corner nodes with the inlet and outlet)
        wall_nodes = []
        for wall_name in ["NoSlip2D_top_wall", "NoSlip2D_bottom_wall"]:
            wall_nodes += list(self.model_part.GetSubModelPart(wall_name).Nodes)
        for node in wall_nodes:
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)

        process = self.__CreateProcess()
        process.ExecuteInitializeSolutionStep()

        self.assertAlmostEqual(self.__CalculateMassFlow(self.inlet_model_part), 2.0, 12)
        self.assertAlmostEqual(self.__CalculateMassFlow(self.outlet_model_part), -2.0, 12)
        for node in wall_nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [0.0, 0.0, 0.0], 12)

        # The walls nodes must remain fixed as these are not fixed by the mass flow process
        process.ExecuteFinalizeSolutionStep()
        wall_node_ids = [node.Id for node in wall_nodes]
        for node in self.model_part.Nodes:
            self.assertEqual(node.IsFixed(KratosMultiphysics.VELOCITY_X), node.Id in wall_node_ids)

    def test_uniform_partially_fixed_nodes(self):
        # Emulate a symmetry condition in the top and bottom walls (only the normal component is fixed)
        # The tangential component of the corner nodes comes from the previous step and must be imposed as well
        wall_nodes = []
        for wall_name in ["NoSlip2D_top_wall", "NoSlip2D_bottom_wall"]:
            wall_nodes += list(self.model_part.GetSubModelPart(wall_name).Nodes)
        for node in wall_nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, [5.0, 0.0, 0.0])
            node.Fix(KratosMultiphysics.VELOCITY_Y)

        process = self.__CreateProcess()
        process.ExecuteInitializeSolutionStep()

        self.assertAlmostEqual(self.__CalculateMassFlow(self.inlet_model_part), 2.0, 12)
        self.assertAlmostEqual(self.__CalculateMassFlow(self.outlet_model_part), -2.0, 12)
        for node in list(self.inlet_model_part.Nodes) + list(self.outlet_model_part.Nodes):
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [2.0 / self.density, 0.0, 0.0], 12)
            self.assertTrue(node.IsFixed(KratosMultiphysics.VELOCITY_X))

        # Only the components fixed by the mass flow process are freed
        process.ExecuteFinalizeSolutionStep()
        wall_node_ids = [node.Id for node in wall_nodes]
        for node in self.model_part.Nodes:
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_X))
            self.assertEqual(node.IsFixed(KratosMultiphysics.VELOCITY_Y), node.Id in wall_node_ids)

    def test_uniform_shared_nodes(self):
        # Split the right wall in two outlets sharing the middle node
        right_wall = self.model_part.GetSubModelPart("NoSlip2D_right_wall")
        for name, is_lower in [("outlet_low", True), ("outlet_high", False)]:
            outlet = self.model_part.CreateSubModelPart(name)
            condition_ids = [condition.Id for condition in right_wall.Conditions if (condition.GetGeometry().Center().Y < 0.6) == is_lower]
            outlet.AddConditions(condition_ids)
            outlet.AddNodes(sorted({node.Id for condition_id in condition_ids for node in self.model_part.GetCondition(condition_id).GetGeometry()}))
        outlet_low = self.model_part.GetSubModelPart("outlet_low")
        outlet_high = self.model_part.GetSubModelPart("outlet_high")
        self.assertEqual(len({node.Id for node in outlet_low.Nodes} & {node.Id for node in outlet_high.Nodes}), 1)

        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "FluidModelPart",
            "boundaries"      : [
                {"model_part_name" : "FluidModelPart.NoSlip2D_left_wall", "mass_flow" :  3.0},
                {"model_part_name" : "FluidModelPart.outlet_low",         "mass_flow" : -1.2},
                {"model_part_name" : "FluidModelPart.outlet_high",        "mass_flow" : -1.8}],
            "density"         : 1.5
        }""")
        process = apply_mass_flow_process.ApplyMassFlowProcess(self.model, settings)
        process.ExecuteInitializeSolutionStep()

        self.assertAlmostEqual(self.__CalculateMassFlow(self.inlet_model_part), 3.0, 12)
        self.assertAlmostEqual(self.__CalculateMassFlow(outlet_low), -1.2, 12)
        self.assertAlmostEqual(self.__CalculateMassFlow(outlet_high), -1.8, 12)

    def test_rescale_existing(self):
        for node in self.inlet_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, [node.Y * (1.0 - node.Y), 0.0, 0.0])

        process = self.__CreateProcess(inlet_profile="rescale_existing")
        process.ExecuteInitializeSolutionStep()

        self.assertAlmostEqual(self.__CalculateMassFlow(self.inlet_model_part), 2.0, 12)
        self.assertAlmostEqual(self.__CalculateMassFlow(self.outlet_model_part), -2.0, 12)

        # Check that the parabolic profile is preserved (0.16 is the unit square flow rate of the interpolated parabola)
        scale_factor = 2.0 / (self.density * 0.16)
        for node in self.inlet_model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), [scale_factor * node.Y * (1.0 - node.Y), 0.0, 0.0], 12)
            self.assertFalse(node.IsFixed(KratosMultiphysics.VELOCITY_X))

    def test_rescale_existing_zero_profile(self):
        process = self.__CreateProcess(inlet_profile="rescale_existing")
        with self.assertRaisesRegex(Exception, "has zero flow rate"):
            process.ExecuteInitializeSolutionStep()

    def test_time_function(self):
        process = self.__CreateProcess(inlet_mass_flow='"2.0*t"', outlet_mass_flow='"-2.0*t"')

        self.model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.5
        process.ExecuteInitializeSolutionStep()
        self.assertAlmostEqual(self.__CalculateMassFlow(self.inlet_model_part), 1.0, 12)
        self.assertAlmostEqual(self.__CalculateMassFlow(self.outlet_model_part), -1.0, 12)
        process.ExecuteFinalizeSolutionStep()

        self.model_part.ProcessInfo[KratosMultiphysics.TIME] = 2.0
        process.ExecuteInitializeSolutionStep()
        self.assertAlmostEqual(self.__CalculateMassFlow(self.inlet_model_part), 4.0, 12)
        self.assertAlmostEqual(self.__CalculateMassFlow(self.outlet_model_part), -4.0, 12)

    def test_density_from_properties(self):
        self.model_part.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY, 4.0)
        self.density = 4.0

        process = self.__CreateProcess(density_settings="")
        process.ExecuteInitializeSolutionStep()

        self.assertAlmostEqual(self.__CalculateMassFlow(self.inlet_model_part), 2.0, 12)
        self.assertAlmostEqual(KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRate(self.inlet_model_part), -0.5, 12)

    def test_check_balance(self):
        process = self.__CreateProcess(outlet_mass_flow="-1.9")
        with self.assertRaisesRegex(Exception, "do not add up to zero"):
            process.ExecuteInitializeSolutionStep()

        process = self.__CreateProcess(outlet_mass_flow="-1.9", check_balance=False)
        process.ExecuteInitializeSolutionStep()
        self.assertAlmostEqual(self.__CalculateMassFlow(self.inlet_model_part), 2.0, 12)
        self.assertAlmostEqual(self.__CalculateMassFlow(self.outlet_model_part), -1.9, 12)

    def __CreateProcess(self, inlet_mass_flow="2.0", outlet_mass_flow="-2.0", inlet_profile="uniform", density_settings='"density" : 1.5,', constrained=True, check_balance=True):
        settings = KratosMultiphysics.Parameters(f"""{{
            "model_part_name" : "FluidModelPart",
            "boundaries"      : [{{
                "model_part_name" : "FluidModelPart.NoSlip2D_left_wall",
                "mass_flow"       : {inlet_mass_flow},
                "profile"         : "{inlet_profile}"
            }},{{
                "model_part_name" : "FluidModelPart.NoSlip2D_right_wall",
                "mass_flow"       : {outlet_mass_flow}
            }}],
            {density_settings}
            "constrained"     : {str(constrained).lower()},
            "check_balance"   : {str(check_balance).lower()}
        }}""")
        return apply_mass_flow_process.ApplyMassFlowProcess(self.model, settings)

    def __CalculateMassFlow(self, boundary_model_part):
        # Note that the flow rate is negative if entering the domain
        return -self.density * KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRate(boundary_model_part)

if __name__ == '__main__':
    KratosUnittest.main()
