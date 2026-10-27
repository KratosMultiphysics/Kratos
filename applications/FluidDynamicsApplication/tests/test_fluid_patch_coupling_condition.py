"""Two independent P1 or Q1 patches coupled across a linear Stokes solution."""

import json
import math

import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as CFD
import KratosMultiphysics.KratosUnittest as KratosUnittest


def exact_velocity(x, y):
    return [1.0 + 2.0 * x + 3.0 * y, -2.0 + 4.0 * x - 2.0 * y, 0.0]


def exact_pressure(x, y):
    return 1.0 + 2.0 * x - 3.0 * y


def create_manufactured_model(coefficient, first_divisions=4,
                              second_divisions=4, element_type="triangle",
                              distorted=False):
    """Create [0,1]x[0,1] and [1,2]x[0,1], sharing no node objects.

    u=(1+2x+3y, -2+4x-2y) is divergence free. For constant viscosity,
    -div(sigma)=grad(p)=(2,-3); density=1 gives BODY_FORCE=(2,-3).
    Only the exterior velocity and one global pressure gauge are fixed.
    """
    element_names = {"triangle": "SymbolicStokes2D3N",
                     "quadrilateral": "SymbolicStokes2D4N"}
    element_name = element_names[element_type]
    model = KM.Model()
    root = model.CreateModelPart("Fluid")
    for variable in (KM.VELOCITY, KM.ACCELERATION, KM.PRESSURE, KM.BODY_FORCE):
        root.AddNodalSolutionStepVariable(variable)
    root.SetBufferSize(3)
    root.ProcessInfo[KM.DOMAIN_SIZE] = 2
    root.ProcessInfo[KM.DELTA_TIME] = 1.0
    root.ProcessInfo[KM.DYNAMIC_TAU] = 0.0
    # SymbolicStokes integrates time internally. Zero BDF coefficients give
    # its steady operator and select the condition's local-system path.
    root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([0.0, 0.0, 0.0])
    next_node_id = next_element_id = 1
    for index, divisions in enumerate((first_divisions, second_divisions)):
        patch = root.CreateSubModelPart(f"Patch_{index}")
        x = float(index)
        domain = KM.Quadrilateral2D4(
            KM.Node(1, x, 0.0, 0.0), KM.Node(2, x, 1.0, 0.0),
            KM.Node(3, x + 1.0, 1.0, 0.0), KM.Node(4, x + 1.0, 0.0, 0.0))
        mesh_settings = KM.Parameters(json.dumps({
            "element_name": element_name,
            "number_of_divisions": [divisions, divisions, 1],
            "create_skin_sub_model_part": False,
            "elements_properties_id": index + 1,
            "start_node_id": next_node_id,
            "start_element_id": next_element_id,
        }))
        KM.StructuredMeshGeneratorProcess(
            domain, patch, mesh_settings).Execute()
        next_node_id += patch.NumberOfNodes()
        next_element_id += patch.NumberOfElements()
        material = root.GetProperties()[index + 1]
        material[KM.DENSITY] = 1.0
        material[KM.DYNAMIC_VISCOSITY] = 1.3
        material[KM.CONSTITUTIVE_LAW] = CFD.Newtonian2DLaw()
        KM.VariableUtils().SetFlag(KM.ACTIVE, True, patch.Elements)

    if distorted:
        # Keep the exterior and interface fixed while making the interior
        # Q1 maps non-affine. Exact linear physical fields remain representable.
        for node in root.Nodes:
            node.X += 0.08 * math.sin(math.pi * node.X) * math.sin(math.pi * node.Y)
            node.X0 = node.X

    for node in root.Nodes:
        for variable in (KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE):
            node.AddDof(variable)
        node.SetSolutionStepValue(KM.BODY_FORCE, [2.0, -3.0, 0.0])
        # The initial interior/interface unknowns are zero, not the exact
        # solution. The interface conditions must recover their values.
        if (abs(node.X) < 1e-12 or abs(node.X - 2.0) < 1e-12 or
                abs(node.Y) < 1e-12 or abs(node.Y - 1.0) < 1e-12):
            node.SetSolutionStepValue(KM.VELOCITY,
                                      exact_velocity(node.X, node.Y))
            node.Fix(KM.VELOCITY_X)
            node.Fix(KM.VELOCITY_Y)
        if abs(node.X) < 1e-12 and abs(node.Y) < 1e-12:
            node.SetSolutionStepValue(KM.PRESSURE,
                                      exact_pressure(node.X, node.Y))
            node.Fix(KM.PRESSURE)

    properties = root.CreateNewProperties(3)
    properties[KM.PENALTY_COEFFICIENT] = coefficient
    settings = KM.Parameters('''{
        "first_model_part_name": "Fluid.Patch_0",
        "second_model_part_name": "Fluid.Patch_1",
        "interface_model_part_name": "Fluid.Coupling",
        "properties_id": 3
    }''')
    CFD.CreateFluidPatchCouplingConditionsProcess(model, settings).Execute()
    return model, root


def solve_manufactured_problem(coefficient, first_divisions=4,
                               second_divisions=4, element_type="triangle",
                               distorted=False):
    """Solve with Kratos elements, conditions, assembly and a direct solver."""
    model, root = create_manufactured_model(
        coefficient, first_divisions, second_divisions, element_type,
        distorted)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    strategy = KM.ResidualBasedLinearStrategy(
        root, scheme, builder, False, False, False, False)
    strategy.SetEchoLevel(0)
    strategy.Initialize()
    strategy.Check()
    strategy.Solve()
    errors = {"velocity": 0.0, "pressure": 0.0, "velocity_jump": 0.0,
              "pressure_jump": 0.0}
    for node in root.Nodes:
        velocity = node.GetSolutionStepValue(KM.VELOCITY)
        exact = exact_velocity(node.X, node.Y)
        errors["velocity"] = max(errors["velocity"],
                                 *(abs(velocity[d] - exact[d])
                                   for d in range(2)))
        errors["pressure"] = max(
            errors["pressure"], abs(node.GetSolutionStepValue(KM.PRESSURE) -
                                    exact_pressure(node.X, node.Y)))

    for condition in root.GetSubModelPart("Coupling").Conditions:
        segment = condition.GetValue(CFD.COUPLING_INTERFACE_COORDINATES)
        y = 0.5 * (segment[0, 1] + segment[1, 1])
        geometry = list(condition.GetGeometry())
        parent_nodes = len(geometry) // 2
        traces = []
        for side in range(2):
            edge = sorted((n for n in geometry[parent_nodes * side:parent_nodes * (side + 1)]
                           if abs(n.X - 1.0) < 1e-12), key=lambda n: n.Y)
            weight = (y - edge[0].Y) / (edge[1].Y - edge[0].Y)
            velocity = [(1.0 - weight) *
                        edge[0].GetSolutionStepValue(KM.VELOCITY)[d] +
                        weight * edge[1].GetSolutionStepValue(KM.VELOCITY)[d]
                        for d in range(2)]
            pressure = ((1.0 - weight) *
                        edge[0].GetSolutionStepValue(KM.PRESSURE) +
                        weight * edge[1].GetSolutionStepValue(KM.PRESSURE))
            traces.append((velocity, pressure))
        errors["velocity_jump"] = max(
            errors["velocity_jump"],
            *(abs(traces[0][0][d] - traces[1][0][d]) for d in range(2)))
        errors["pressure_jump"] = max(
            errors["pressure_jump"], abs(traces[0][1] - traces[1][1]))
    strategy.Clear()
    return model, root, errors


class FluidPatchCouplingConditionTest(KratosUnittest.TestCase):

    def _check_manufactured(self, coefficient, divisions,
                            element_type="triangle", distorted=False):
        model, root, errors = solve_manufactured_problem(
            coefficient, *divisions, element_type=element_type,
            distorted=distorted)
        for name, error in errors.items():
            self.assertTrue(math.isfinite(error), name)
            self.assertLess(error, 1e-8, f"{name}: {errors}")
        first = root.GetSubModelPart("Patch_0")
        second = root.GetSubModelPart("Patch_1")
        self.assertTrue({n.Id for n in first.Nodes}.isdisjoint(
            {n.Id for n in second.Nodes}))
        self.assertEqual(sum(n.IsFixed(KM.PRESSURE) for n in root.Nodes), 1)
        interface = root.GetSubModelPart("Coupling")
        expected_segments = sum(divisions) - math.gcd(*divisions)
        self.assertEqual(interface.NumberOfConditions(), expected_segments)
        total_length = 0.0
        free_interface_nodes = []
        for node in root.Nodes:
            if abs(node.X - 1.0) < 1e-12 and 1e-12 < node.Y < 1.0 - 1e-12:
                self.assertFalse(node.IsFixed(KM.VELOCITY_X))
                self.assertFalse(node.IsFixed(KM.VELOCITY_Y))
                self.assertFalse(node.IsFixed(KM.PRESSURE))
                free_interface_nodes.append(node.Id)
        self.assertGreater(len(free_interface_nodes), 0)
        parent_nodes = 3 if element_type == "triangle" else 4
        for condition in interface.Conditions:
            self.assertEqual(condition.GetGeometry().PointsNumber(),
                             2 * parent_nodes)
            segment = condition.GetValue(
                CFD.COUPLING_INTERFACE_COORDINATES)
            self.assertAlmostEqual(segment[0, 0], 1.0)
            self.assertAlmostEqual(segment[1, 0], 1.0)
            total_length += segment[1, 1] - segment[0, 1]
        self.assertAlmostEqual(total_length, 1.0)
        KM.Logger.PrintInfo("FluidPatchCouplingManufactured",
                            f"coefficient={coefficient}, "
                            f"divisions={divisions}, element_type={element_type}, "
                            f"distorted={distorted}, "
                            f"errors={errors}")

    def test_linear_manufactured_penalty(self):
        self._check_manufactured(40.0, (4, 4))

    def test_linear_manufactured_penalty_free(self):
        self._check_manufactured(0.0, (4, 4))

    def test_nonmatching_linear_manufactured_penalty(self):
        self._check_manufactured(40.0, (4, 6))

    def test_nonmatching_linear_manufactured_penalty_free(self):
        self._check_manufactured(0.0, (4, 6))

    def test_quadrilateral_linear_manufactured_penalty(self):
        self._check_manufactured(40.0, (4, 4), "quadrilateral")

    def test_quadrilateral_linear_manufactured_penalty_free(self):
        self._check_manufactured(0.0, (4, 4), "quadrilateral")

    def test_nonmatching_quadrilateral_linear_manufactured_penalty(self):
        self._check_manufactured(40.0, (4, 6), "quadrilateral")

    def test_nonmatching_quadrilateral_linear_manufactured_penalty_free(self):
        self._check_manufactured(0.0, (4, 6), "quadrilateral")

    def test_distorted_quadrilateral_linear_manufactured(self):
        for coefficient in (40.0, 0.0):
            with self.subTest(coefficient=coefficient):
                self._check_manufactured(coefficient, (4, 6),
                                         "quadrilateral", distorted=True)

    def test_linear_field_balances_volume_and_interface_residuals(self):
        self._check_volume_and_interface_residuals("triangle")

    def test_quadrilateral_field_balances_volume_and_interface_residuals(self):
        self._check_volume_and_interface_residuals("quadrilateral")

    def _check_volume_and_interface_residuals(self, element_type):
        _, root = create_manufactured_model(0.0, 4, 6, element_type)
        for node in root.Nodes:
            node.SetSolutionStepValue(KM.VELOCITY,
                                      exact_velocity(node.X, node.Y))
            node.SetSolutionStepValue(KM.PRESSURE,
                                      exact_pressure(node.X, node.Y))
        residual = {n.Id: [0.0, 0.0, 0.0] for n in root.Nodes}

        def add_residual(entity):
            entity.Initialize(root.ProcessInfo)
            lhs, rhs = KM.Matrix(), KM.Vector()
            entity.CalculateLocalSystem(lhs, rhs, root.ProcessInfo)
            for i, node in enumerate(entity.GetGeometry()):
                for d in range(3):
                    residual[node.Id][d] += rhs[3 * i + d]

        for element in root.Elements:
            add_residual(element)
        self.assertGreater(max(abs(residual[n.Id][d]) for n in root.Nodes
                               if abs(n.X - 1.0) < 1e-12 and
                               1e-12 < n.Y < 1.0 - 1e-12
                               for d in range(2)), 0.1)
        for condition in root.GetSubModelPart("Coupling").Conditions:
            add_residual(condition)
        for node in root.Nodes:
            for d, variable in enumerate((KM.VELOCITY_X, KM.VELOCITY_Y,
                                           KM.PRESSURE)):
                if not node.IsFixed(variable):
                    self.assertAlmostEqual(residual[node.Id][d], 0.0,
                                           delta=1e-10)

    def test_creator_rejects_duplicate_execution(self):
        model, root = create_manufactured_model(0.0)
        settings = KM.Parameters('''{
            "first_model_part_name": "Fluid.Patch_0",
            "second_model_part_name": "Fluid.Patch_1",
            "interface_model_part_name": "Fluid.Coupling",
            "properties_id": 3
        }''')
        with self.assertRaisesRegex(RuntimeError,
                                    "execute the creator only once"):
            CFD.CreateFluidPatchCouplingConditionsProcess(
                model, settings).Execute()


if __name__ == "__main__":
    KratosUnittest.main()
