"""Initial-state stress assignment for the supported 2D and 3D Dam geometries.

Source and destination always share an identical mesh; only the Gauss-point
stresses are extracted and assigned through the lightweight initial-state
utility used by the JSON/gzip import/export processes.
"""
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.DamApplication as Dam
import KratosMultiphysics.KratosUnittest as KratosUnittest


GEOMETRIES = {
    'T3': {
        'element': 'SmallDisplacementSolidElement2D3N', 'dim': 2, 'gps': 1,
        'coords': [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)]},
    'T6': {
        'element': 'SmallDisplacementSolidElement2D6N', 'dim': 2, 'gps': 3,
        'coords': [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
                   (0.5, 0., 0.), (0.5, 0.5, 0.), (0., 0.5, 0.)]},
    'Q4': {
        'element': 'SmallDisplacementSolidElement2D4N', 'dim': 2, 'gps': 4,
        'coords': [(0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.)]},
    'Q8': {
        'element': 'SmallDisplacementSolidElement2D8N', 'dim': 2, 'gps': 9,
        'coords': [(0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
                   (0.5, 0., 0.), (1., 0.5, 0.), (0.5, 1., 0.), (0., 0.5, 0.)]},
    'Q9': {
        'element': 'SmallDisplacementSolidElement2D9N', 'dim': 2, 'gps': 9,
        'coords': [(0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
                   (0.5, 0., 0.), (1., 0.5, 0.), (0.5, 1., 0.), (0., 0.5, 0.),
                   (0.5, 0.5, 0.)]},
    'T4': {
        'element': 'SmallDisplacementSolidElement3D4N', 'dim': 3, 'gps': 1,
        'coords': [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.)]},
    'T10': {
        'element': 'SmallDisplacementSolidElement3D10N', 'dim': 3, 'gps': 4,
        'coords': [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.),
                   (0.5, 0., 0.), (0.5, 0.5, 0.), (0., 0.5, 0.),
                   (0., 0., 0.5), (0.5, 0., 0.5), (0., 0.5, 0.5)]},
    'P6': {
        'element': 'SmallDisplacementSolidElement3D6N', 'dim': 3, 'gps': 6,
        'coords': [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
                   (0., 0., 1.), (1., 0., 1.), (0., 1., 1.)]},
    'P15': {
        'element': 'SmallDisplacementSolidElement3D15N', 'dim': 3, 'gps': 9,
        'coords': [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
                   (0., 0., 1.), (1., 0., 1.), (0., 1., 1.),
                   (0.5, 0., 0.), (0.5, 0.5, 0.), (0., 0.5, 0.),
                   (0., 0., 0.5), (1., 0., 0.5), (0., 1., 0.5),
                   (0.5, 0., 1.), (0.5, 0.5, 1.), (0., 0.5, 1.)]},
    'H8': {
        'element': 'SmallDisplacementSolidElement3D8N', 'dim': 3, 'gps': 8,
        'coords': [(0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
                   (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.)]},
    'H20': {
        'element': 'SmallDisplacementSolidElement3D20N', 'dim': 3, 'gps': 27,
        'coords': [(0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
                   (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
                   (0.5, 0., 0.), (1., 0.5, 0.), (0.5, 1., 0.), (0., 0.5, 0.),
                   (0., 0., 0.5), (1., 0., 0.5), (1., 1., 0.5), (0., 1., 0.5),
                   (0.5, 0., 1.), (1., 0.5, 1.), (0.5, 1., 1.), (0., 0.5, 1.)]},
    'H27': {
        'element': 'SmallDisplacementSolidElement3D27N', 'dim': 3, 'gps': 27,
        'coords': [(0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
                   (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
                   (0.5, 0., 0.), (1., 0.5, 0.), (0.5, 1., 0.), (0., 0.5, 0.),
                   (0., 0., 0.5), (1., 0., 0.5), (1., 1., 0.5), (0., 1., 0.5),
                   (0.5, 0., 1.), (1., 0.5, 1.), (0.5, 1., 1.), (0., 0.5, 1.),
                   (0.5, 0.5, 0.), (0.5, 0., 0.5), (1., 0.5, 0.5), (0.5, 1., 0.5),
                   (0., 0.5, 0.5), (0.5, 0.5, 1.), (0.5, 0.5, 0.5)]},
}

# Used only to check the unsupported-geometry rejection path.
UNSUPPORTED_T10_2D = {
    'element': 'SmallDisplacementElement2D10N', 'dim': 2, 'gps': 6,
    'coords': [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
               (1. / 3., 0., 0.), (2. / 3., 0., 0.),
               (2. / 3., 1. / 3., 0.), (1. / 3., 2. / 3., 0.),
               (0., 2. / 3., 0.), (0., 1. / 3., 0.),
               (1. / 3., 1. / 3., 0.)]}

DISPLACEMENT_DOFS = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
REACTION_DOFS = (KM.REACTION_X, KM.REACTION_Y, KM.REACTION_Z)


class GeometryFixture:

    @staticmethod
    def _law(dim, thermal=False, plane_stress=False):
        if dim == 2:
            if thermal:
                return (Dam.ThermalLinearElastic2DPlaneStress() if plane_stress
                        else Dam.ThermalLinearElastic2DPlaneStrain())
            return (SMA.LinearElasticPlaneStress2DLaw() if plane_stress
                    else SMA.LinearElasticPlaneStrain2DLaw())
        return Dam.ThermalLinearElastic3DLaw() if thermal else SMA.LinearElastic3DLaw()

    @staticmethod
    def _build(geometry, law, part_name='MainModelPart', element_id=1,
               connectivity=None, integration_order=None, body=None, fixed_ids=()):
        dim = geometry['dim']
        coords = geometry['coords']
        model = KM.Model()
        part = model.CreateModelPart(part_name, 2)
        for variable in (KM.DISPLACEMENT, KM.REACTION, KM.VOLUME_ACCELERATION,
                         KM.CAUCHY_STRESS_VECTOR, KM.GREEN_LAGRANGE_STRAIN_VECTOR,
                         KM.TEMPERATURE, Dam.NODAL_REFERENCE_TEMPERATURE):
            part.AddNodalSolutionStepVariable(variable)
        part.ProcessInfo[KM.DOMAIN_SIZE] = dim
        part.ProcessInfo[KM.DELTA_TIME] = 1.0
        part.ProcessInfo[KM.TIME] = 0.0

        nodes = []
        for node_id, (x, y, z) in enumerate(coords, 1):
            node = part.CreateNewNode(node_id, x, y, z)
            for dof, reaction in zip(DISPLACEMENT_DOFS, REACTION_DOFS):
                node.AddDof(dof, reaction)
            node.SetSolutionStepValue(KM.TEMPERATURE, 20.0)
            node.SetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE, 20.0)
            if body is not None:
                node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, body)
            if node_id in fixed_ids:
                for dof in DISPLACEMENT_DOFS:
                    node.Fix(dof)
            nodes.append(node)

        properties = part.CreateNewProperties(1)
        properties[KM.YOUNG_MODULUS] = 1.e4
        properties[KM.POISSON_RATIO] = 0.2
        properties[KM.DENSITY] = 2.0
        properties[Dam.THERMAL_EXPANSION] = 1.e-5
        if dim == 2:
            properties[KM.THICKNESS] = 1.0
        if integration_order is not None:
            properties[KM.KratosGlobals.GetVariable('INTEGRATION_ORDER')] = integration_order
        properties[KM.CONSTITUTIVE_LAW] = law

        if connectivity is None:
            connectivity = list(range(1, len(coords) + 1))
        element = part.CreateNewElement(
            geometry['element'], element_id, connectivity, properties)
        element.Initialize(part.ProcessInfo)
        return model, part, nodes, element

    @staticmethod
    def _set_quadratic_displacement(nodes, coords, dim):
        for node, (x, y, z) in zip(nodes, coords):
            if dim == 2:
                node.SetSolutionStepValue(KM.DISPLACEMENT, [
                    0.01 * x + 0.004 * y + 0.002 * x * x + 0.001 * x * y,
                    -0.003 * x + 0.006 * y + 0.0015 * y * y + 0.0005 * x * y,
                    0.0])
            else:
                node.SetSolutionStepValue(KM.DISPLACEMENT, [
                    0.01 * x + 0.004 * y + 0.002 * z +
                    0.001 * x * x + 0.0005 * x * y + 0.0005 * x * z,
                    -0.003 * x + 0.006 * y + 0.001 * z +
                    0.001 * y * y + 0.0005 * x * y + 0.0005 * y * z,
                    0.002 * x - 0.001 * y + 0.005 * z +
                    0.001 * z * z + 0.0005 * x * z + 0.0005 * y * z])

    @staticmethod
    def _stresses(element, part):
        return [list(value) for value in element.CalculateOnIntegrationPoints(
            KM.CAUCHY_STRESS_VECTOR, part.ProcessInfo)]

    def _assert_close(self, actual, expected, tolerance=1.e-10):
        self.assertLessEqual(abs(actual - expected),
                             tolerance * (1.0 + abs(expected)))

    def _assert_stresses_close(self, actual, expected, tolerance=1.e-10):
        self.assertEqual(len(actual), len(expected))
        for actual_point, expected_point in zip(actual, expected):
            self.assertEqual(len(actual_point), len(expected_point))
            for a, e in zip(actual_point, expected_point):
                self._assert_close(a, e, tolerance)

    def _assert_zero_stresses(self, element, part):
        for point in self._stresses(element, part):
            for value in point:
                self.assertEqual(value, 0.0)

    @staticmethod
    def _assign(source, destination):
        stresses = Dam.DamInitialStateTransferUtility.GetFinalCauchyStresses(source)
        Dam.DamInitialStateTransferUtility.AssignInitialStresses(
            destination, stresses)

    @staticmethod
    def _make_strategy(part):
        builder = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        strategy = KM.ResidualBasedLinearStrategy(
            part, scheme, builder, True, False, False, False)
        strategy.SetEchoLevel(0)
        strategy.Initialize()
        strategy.Check()
        return strategy

    def _free_residual(self, part, element, dim):
        residual = [0.0] * (dim * part.NumberOfNodes())
        left_hand_side, right_hand_side = KM.Matrix(), KM.Vector()
        element.CalculateLocalSystem(left_hand_side, right_hand_side, part.ProcessInfo)
        for i, node in enumerate(element.GetGeometry()):
            for d in range(dim):
                residual[dim * (node.Id - 1) + d] += right_hand_side[dim * i + d]
        return residual

class TestDamInitialStateGeometry(GeometryFixture, KratosUnittest.TestCase):

    def test_stress_assignment_all_geometries(self):
        for name, geometry in GEOMETRIES.items():
            with self.subTest(geometry=name):
                dim = geometry['dim']
                law = self._law(dim)
                source_model, source, source_nodes, source_element = self._build(
                    geometry, law)
                self._set_quadratic_displacement(
                    source_nodes, geometry['coords'], dim)
                source_stress = self._stresses(source_element, source)
                self.assertEqual(len(source_stress), geometry['gps'])
                self.assertEqual(len(source_stress[0]), 3 if dim == 2 else 6)

                if dim == 2:
                    self.assertGreater(
                        max(abs(point[2]) for point in source_stress), 1.e-6)
                else:
                    self.assertGreater(max(
                        abs(point[3]) + abs(point[4]) + abs(point[5])
                        for point in source_stress), 1.e-6)
                if geometry['gps'] > 1:
                    self.assertGreater(max(
                        abs(point[j] - source_stress[0][j])
                        for point in source_stress
                        for j in range(len(point))), 1.e-6)

                destination_model, destination, destination_nodes, destination_element = \
                    self._build(geometry, law)
                for node in destination_nodes:
                    for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                        self.assertEqual(value, 0.0)

                self._assign(source, destination)
                destination_stress = self._stresses(
                    destination_element, destination)
                self._assert_stresses_close(destination_stress, source_stress)
                for node in destination_nodes:
                    for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                        self.assertEqual(value, 0.0)
                    self.assertEqual(
                        node.GetSolutionStepValue(KM.TEMPERATURE), 20.0)
                    self.assertEqual(node.GetSolutionStepValue(
                        Dam.NODAL_REFERENCE_TEMPERATURE), 20.0)

    def test_plane_stress_and_thermal_laws(self):
        cases = (
            ('Q4_plane_stress', 'Q4', False, False),
            ('Q4_thermal_plane_strain', 'Q4', True, False),
            ('Q4_thermal_plane_stress', 'Q4', True, True),
            ('T6_thermal_plane_strain', 'T6', True, False),
            ('Q8_thermal_plane_strain', 'Q8', True, False),
            ('Q9_thermal_plane_strain', 'Q9', True, False),
            ('H8_thermal_3d', 'H8', True, False),
            ('T10_thermal_3d', 'T10', True, False),
            ('H20_thermal_3d', 'H20', True, False),
            ('H27_thermal_3d', 'H27', True, False),
        )
        for label, geometry_name, thermal, plane_stress in cases:
            with self.subTest(case=label):
                geometry = GEOMETRIES[geometry_name]
                dim = geometry['dim']
                law = self._law(dim, thermal=thermal, plane_stress=plane_stress)
                source_model, source, source_nodes, source_element = self._build(
                    geometry, law)
                self._set_quadratic_displacement(
                    source_nodes, geometry['coords'], dim)
                source_stress = self._stresses(source_element, source)

                destination_model, destination, _, destination_element = self._build(
                    geometry, law)
                thermal_reference_model, thermal_reference, _, thermal_reference_element = \
                    self._build(geometry, law)
                self._assign(source, destination)
                self._assert_stresses_close(
                    self._stresses(destination_element, destination), source_stress)

                for node in destination.Nodes:
                    node.SetSolutionStepValue(KM.TEMPERATURE, 30.0)
                for node in thermal_reference.Nodes:
                    node.SetSolutionStepValue(KM.TEMPERATURE, 30.0)
                reference_stress = self._stresses(
                    thermal_reference_element, thermal_reference)
                expected_hot = [[a + b for a, b in zip(source_point, reference_point)]
                                for source_point, reference_point in
                                zip(source_stress, reference_stress)]
                self._assert_stresses_close(
                    self._stresses(destination_element, destination), expected_hot)
                for node in destination.Nodes:
                    self.assertEqual(node.GetSolutionStepValue(
                        Dam.NODAL_REFERENCE_TEMPERATURE), 20.0)

    def test_non_default_integration_rules(self):
        cases = (
            ('Q4_order3', 'Q4', 3, 9),
            ('H8_order3', 'H8', 3, 27),
        )
        for label, geometry_name, order, expected_points in cases:
            with self.subTest(case=label):
                geometry = GEOMETRIES[geometry_name]
                dim = geometry['dim']
                law = self._law(dim)
                source_model, source, source_nodes, source_element = self._build(
                    geometry, law, integration_order=order)
                self._set_quadratic_displacement(
                    source_nodes, geometry['coords'], dim)
                source_stress = self._stresses(source_element, source)
                self.assertEqual(len(source_stress), expected_points)
                self.assertGreater(max(
                    abs(point[j] - source_stress[0][j])
                    for point in source_stress
                    for j in range(len(point))), 1.e-6)

                destination_model, destination, destination_nodes, destination_element = \
                    self._build(geometry, law, integration_order=order)
                self._assign(source, destination)
                self._assert_stresses_close(
                    self._stresses(destination_element, destination), source_stress)
                for node in destination_nodes:
                    for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                        self.assertEqual(value, 0.0)

    def test_source_mesh_motion_stress_extraction(self):
        geometry = GEOMETRIES['Q4']
        law = self._law(2)
        source_model, source, source_nodes, source_element = self._build(
            geometry, law)
        self._set_quadratic_displacement(source_nodes, geometry['coords'], 2)
        for node in source_nodes:
            displacement = node.GetSolutionStepValue(KM.DISPLACEMENT)
            node.X += displacement[0]
            node.Y += displacement[1]
            node.Z += displacement[2]
        source_stress = self._stresses(source_element, source)

        destination_model, destination, _, destination_element = self._build(
            geometry, law)
        self._assign(source, destination)
        self._assert_stresses_close(
            self._stresses(destination_element, destination), source_stress)

    def test_equilibrium_and_incremental_response_all_geometries(self):
        fixed_ids_2d = (1, 2)
        fixed_ids_3d = (1, 2, 3)
        for name, geometry in GEOMETRIES.items():
            with self.subTest(geometry=name):
                dim = geometry['dim']
                body = [0.0, -10.0, 0.0] if dim == 2 else [0.0, 0.0, -10.0]
                fixed_ids = fixed_ids_2d if dim == 2 else (
                    (1, 2, 3, 4) if name in ('H8', 'H20', 'H27') else fixed_ids_3d)
                law = self._law(dim)

                source_model, source, _, source_element = self._build(
                    geometry, law, body=body, fixed_ids=fixed_ids)
                source_strategy = self._make_strategy(source)
                source_strategy.Solve()
                gravity_stress = self._stresses(source_element, source)
                gravity_displacement = [
                    list(node.GetSolutionStepValue(KM.DISPLACEMENT))
                    for node in source.Nodes]
                gravity_reaction = [
                    list(node.GetSolutionStepValue(KM.REACTION))
                    for node in source.Nodes]
                self.assertGreater(max(abs(value) for displacement in gravity_displacement
                                       for value in displacement), 1.e-8)
                self.assertGreater(max(abs(value) for point in gravity_stress
                                       for value in point), 1.e-6)

                destination_model, destination, destination_nodes, destination_element = \
                    self._build(geometry, law, body=body, fixed_ids=fixed_ids)
                free_dofs = sum(
                    1 for node in destination.Nodes for d in range(dim)
                    if not node.IsFixed(DISPLACEMENT_DOFS[d]))
                self.assertGreater(free_dofs, 0)
                destination_strategy = self._make_strategy(destination)
                self._assign(source, destination)
                destination_strategy.Solve()
                self.assertLess(max(abs(value) for node in destination_nodes
                                    for value in node.GetSolutionStepValue(KM.DISPLACEMENT)),
                                1.e-12)
                self._assert_stresses_close(
                    self._stresses(destination_element, destination), gravity_stress)
                destination_reaction = [
                    list(node.GetSolutionStepValue(KM.REACTION))
                    for node in destination.Nodes]
                for expected, actual in zip(gravity_reaction, destination_reaction):
                    for e, a in zip(expected, actual):
                        self._assert_close(a, e)

                residual = self._free_residual(destination, destination_element, dim)
                for node in destination.Nodes:
                    for d in range(dim):
                        if not node.IsFixed(DISPLACEMENT_DOFS[d]):
                            self.assertLess(abs(residual[dim * (node.Id - 1) + d]), 1.e-9)

                measure = (destination_element.GetGeometry().Area() if dim == 2
                           else destination_element.GetGeometry().Volume())
                if dim == 2:
                    expected_reaction = -2.0 * measure * 1.0 * body[1]
                    actual_reaction = sum(node.GetSolutionStepValue(KM.REACTION)[1]
                                          for node in destination.Nodes)
                else:
                    expected_reaction = -2.0 * measure * body[2]
                    actual_reaction = sum(node.GetSolutionStepValue(KM.REACTION)[2]
                                          for node in destination.Nodes)
                self._assert_close(actual_reaction, expected_reaction, 1.e-9)

                increment = [2.0, -10.0, 0.0] if dim == 2 else [2.0, 0.0, -10.0]
                for node in source.Nodes:
                    node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, increment)
                for node in destination.Nodes:
                    node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, increment)
                source_strategy.Solve()
                destination_strategy.Solve()
                source_displacement = [
                    list(node.GetSolutionStepValue(KM.DISPLACEMENT))
                    for node in source.Nodes]
                source_increment = [
                    a - b for after, before in
                    zip(source_displacement, gravity_displacement)
                    for a, b in zip(after, before)]
                destination_increment = [
                    value for node in destination.Nodes
                    for value in node.GetSolutionStepValue(KM.DISPLACEMENT)]
                self.assertGreater(max(abs(value) for value in source_increment), 1.e-8)
                self.assertGreater(max(abs(value) for value in destination_increment), 1.e-8)
                increment_stress = self._stresses(source_element, source)
                self.assertGreater(max(
                    abs(point[j] - gravity_stress[i][j])
                    for i, point in enumerate(increment_stress)
                    for j in range(len(point))), 1.e-6)
                for expected, actual in zip(source_increment, destination_increment):
                    self._assert_close(actual, expected)
                self._assert_stresses_close(
                    self._stresses(destination_element, destination),
                    self._stresses(source_element, source))
                for expected, actual in zip(
                        [node.GetSolutionStepValue(KM.REACTION) for node in source.Nodes],
                        [node.GetSolutionStepValue(KM.REACTION)
                         for node in destination.Nodes]):
                    for e, a in zip(expected, actual):
                        self._assert_close(a, e)

    def test_rejects_unsupported_assignment_configurations(self):
        geometry = GEOMETRIES['Q4']
        law = self._law(2)
        source_model, source, _, _ = self._build(geometry, law)

        destination_model, destination, _, destination_element = self._build(
            geometry, Dam.ThermalSimoJuLocalDamagePlaneStrain2DLaw())
        with self.assertRaisesRegex(RuntimeError, 'unsupported constitutive law'):
            self._assign(source, destination)

        destination_model, destination, _, destination_element = self._build(
            geometry, law)
        destination.Nodes[1].SetSolutionStepValue(
            KM.DISPLACEMENT, [0.1, 0.0, 0.0])
        with self.assertRaisesRegex(RuntimeError, 'nonzero displacement'):
            self._assign(source, destination)

        destination_model, destination, _, destination_element = self._build(
            geometry, law)
        self._assign(source, destination)
        with self.assertRaisesRegex(RuntimeError, 'already has an initial state'):
            self._assign(source, destination)

        destination_model, destination, _, destination_element = self._build(
            geometry, law)
        destination.CreateNewElement(
            geometry['element'], 2, [1, 2, 3, 4], destination.Properties[1])
        with self.assertRaisesRegex(RuntimeError, 'imported stress map contains'):
            self._assign(source, destination)

        unsupported_source_model, unsupported_source, _, _ = self._build(
            UNSUPPORTED_T10_2D, self._law(2))
        unsupported_destination_model, unsupported_destination, _, unsupported_element = \
            self._build(UNSUPPORTED_T10_2D, self._law(2))
        with self.assertRaisesRegex(RuntimeError, 'unsupported geometry'):
            self._assign(unsupported_source, unsupported_destination)
        self._assert_zero_stresses(unsupported_element, unsupported_destination)

if __name__ == '__main__':
    KratosUnittest.main()
