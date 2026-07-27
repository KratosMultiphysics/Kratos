import math

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestPatchTestCuttingPattern(KratosUnittest.TestCase):

    def _add_variables(self, mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

    def _add_dofs(self, mp):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

    def _apply_material_properties(self, mp):
        props = mp.GetProperties()[1]
        props.SetValue(KratosMultiphysics.THICKNESS, 0.01)
        props.SetValue(KratosMultiphysics.YOUNG_MODULUS, 1.0e6)
        props.SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
        cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
        props.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, cl)

    def _create_shallow_pyramid_patch(self, mp):
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 1.0, 1.0, 0.0)
        mp.CreateNewNode(4, 0.0, 1.0, 0.0)
        mp.CreateNewNode(5, 0.5, 0.5, 0.2)

        element_name = "MembraneCuttingPatternElement3D3N"
        mp.CreateNewElement(element_name, 1, [1, 2, 5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [2, 3, 5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 3, [3, 4, 5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [4, 1, 5], mp.GetProperties()[1])

    def _apply_BCs(self, mp):
        # Keep the whole solve in-plane (this is a planar flattening problem)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)
        mp.Nodes[1].Fix(KratosMultiphysics.DISPLACEMENT_X)
        mp.Nodes[1].Fix(KratosMultiphysics.DISPLACEMENT_Y)
        mp.Nodes[2].Fix(KratosMultiphysics.DISPLACEMENT_Y)

    def _setup_cutting_pattern_strategy(self, mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = StructuralMechanicsApplication.CuttingPatternScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-10, 1e-10)
        convergence_criterion.SetEchoLevel(0)

        initial_flattening_settings = KratosMultiphysics.Parameters("""{
            "projection_type"  : "planar_fixed_direction",
            "global_direction" : [0.0, 0.0, 1.0],
            "echo_level"       : 0
        }""")

        max_iters = 100
        compute_reactions = True
        reform_step_dofs = True
        move_mesh_flag = False
        use_relaxation = False
        strategy = StructuralMechanicsApplication.CuttingPatternStrategy(
            mp, scheme, convergence_criterion, builder_and_solver,
            mp, False, "none", initial_flattening_settings,
            max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, use_relaxation)
        strategy.SetEchoLevel(0)
        return strategy

    def test_cutting_pattern_shallow_pyramid(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("cutting_pattern_part")
        mp.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_shallow_pyramid_patch(mp)
        self._add_dofs(mp)
        self._apply_BCs(mp)

        # snapshot reference geometry before solving
        ref_coords = {node.Id: (node.X0, node.Y0, node.Z0) for node in mp.Nodes}

        strategy = self._setup_cutting_pattern_strategy(mp)
        strategy.Initialize()
        strategy.Check()
        strategy.Solve()

        for node in mp.Nodes:
            ref_x0, ref_y0, ref_z0 = ref_coords[node.Id]
            self.assertAlmostEqual(node.X0, ref_x0)
            self.assertAlmostEqual(node.Y0, ref_y0)
            self.assertAlmostEqual(node.Z0, ref_z0)

        # the flattened shape must land exactly in the global Z=0 plane
        for node in mp.Nodes:
            self.assertAlmostEqual(node.Z0 + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)

        def in_plane_position(node):
            return (
                node.X0 + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),
                node.Y0 + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),
            )

        def distance(id_a, id_b):
            xa, ya = in_plane_position(mp.Nodes[id_a])
            xb, yb = in_plane_position(mp.Nodes[id_b])
            return math.sqrt((xa - xb) ** 2 + (ya - yb) ** 2)

        for id_a, id_b in [(1, 2), (2, 3), (3, 4), (4, 1)]:
            self.assertAlmostEqual(distance(id_a, id_b), 1.0, delta=0.05)

        for corner_id in [1, 2, 3, 4]:
            radial_length = distance(corner_id, 5)
            self.assertGreater(radial_length, 0.70)
            self.assertLess(radial_length, 0.75)

    def test_initial_flattening_utility_fixed_direction(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("initial_flattening_part")
        mp.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_shallow_pyramid_patch(mp)
        self._add_dofs(mp)

        settings = KratosMultiphysics.Parameters("""{
            "projection_type"  : "planar_fixed_direction",
            "global_direction" : [0.0, 0.0, 1.0],
            "echo_level"       : 0
        }""")
        StructuralMechanicsApplication.InitialFlatteningUtility.Execute(mp, settings)

        for node_id, expected_z in [(1, 0.0), (2, 0.0), (3, 0.0), (4, 0.0), (5, -0.2)]:
            node = mp.Nodes[node_id]
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), 0.0)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), expected_z)

    def test_initial_flattening_utility_per_part_mean_normal(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("per_part_flattening")
        mp.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self._add_variables(mp)
        self._apply_material_properties(mp)

        node_coords = {
            1: (0.0, 0.0, 0.0),
            2: (1.0, 0.0, 0.0),
            3: (0.0, 1.0, 1.0),
            4: (2.0, 0.0, 0.0),
            5: (1.0, 1.0, -1.0),
        }
        for nid, (x, y, z) in node_coords.items():
            mp.CreateNewNode(nid, x, y, z)
        self._add_dofs(mp)

        element_name = "MembraneCuttingPatternElement3D3N"
        elem_a = mp.CreateNewElement(element_name, 1, [1, 2, 3], mp.GetProperties()[1])
        elem_b = mp.CreateNewElement(element_name, 2, [2, 4, 5], mp.GetProperties()[1])

        part_a = mp.CreateSubModelPart("PartA")
        part_a.AddNodes([1, 2, 3])
        part_a.AddElements([1])

        part_b = mp.CreateSubModelPart("PartB")
        part_b.AddNodes([2, 4, 5])
        part_b.AddElements([2])

        settings = KratosMultiphysics.Parameters("""{
            "projection_type" : "planar_mean_normal",
            "sub_model_parts" : ["PartA", "PartB"],
            "echo_level"      : 0
        }""")
        StructuralMechanicsApplication.InitialFlatteningUtility.Execute(mp, settings)

        def unit_normal(p0, p1, p2):
            v1 = [p1[i] - p0[i] for i in range(3)]
            v2 = [p2[i] - p0[i] for i in range(3)]
            cross = [
                v1[1] * v2[2] - v1[2] * v2[1],
                v1[2] * v2[0] - v1[0] * v2[2],
                v1[0] * v2[1] - v1[1] * v2[0],
            ]
            norm = math.sqrt(sum(c * c for c in cross))
            return [c / norm for c in cross]

        normal_a = unit_normal(node_coords[1], node_coords[2], node_coords[3])
        normal_b = unit_normal(node_coords[2], node_coords[4], node_coords[5])
        shared_sum = [normal_a[i] + normal_b[i] for i in range(3)]
        shared_norm = math.sqrt(sum(c * c for c in shared_sum))
        normal_shared = [c / shared_norm for c in shared_sum]

        def expected_proj(pos, normal):
            dot = sum(pos[i] * normal[i] for i in range(3))
            return [-dot * normal[i] for i in range(3)]

        expected = {
            1: expected_proj(node_coords[1], normal_a),
            2: expected_proj(node_coords[2], normal_shared),
            3: expected_proj(node_coords[3], normal_a),
            4: expected_proj(node_coords[4], normal_b),
            5: expected_proj(node_coords[5], normal_b),
        }

        for node_id, expected_disp in expected.items():
            disp = mp.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            for i in range(3):
                self.assertAlmostEqual(disp[i], expected_disp[i], places=10)


if __name__ == '__main__':
    KratosUnittest.main()
