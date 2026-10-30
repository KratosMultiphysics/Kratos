import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as Dam
import KratosMultiphysics.PoromechanicsApplication as Poro
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestDamElasticInitialStateTransfer(KratosUnittest.TestCase):
    """Two-stage gravity and horizontal body-force loading on a 1x2 body."""

    def _body(self, quadrilaterals=False, thermal=False):
        model = KM.Model()
        part = model.CreateModelPart("body", 2)
        for variable in (KM.DISPLACEMENT, KM.REACTION, KM.VOLUME_ACCELERATION,
                         KM.TEMPERATURE, Dam.NODAL_REFERENCE_TEMPERATURE,
                         Poro.NODAL_CAUCHY_STRESS_TENSOR,
                         Dam.INITIAL_NODAL_CAUCHY_STRESS_TENSOR):
            part.AddNodalSolutionStepVariable(variable)
        part.ProcessInfo[KM.DOMAIN_SIZE] = 2
        part.ProcessInfo[KM.DELTA_TIME] = 1.0
        for node_id, x, y in ((1, 0., 0.), (2, 1., 0.), (3, 0., 1.),
                              (4, 1., 1.), (5, 0., 2.), (6, 1., 2.)):
            node = part.CreateNewNode(node_id, x, y, 0.)
            for dof, reaction in ((KM.DISPLACEMENT_X, KM.REACTION_X),
                                  (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                                  (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
                node.AddDof(dof, reaction)
            if y == 0.:
                node.Fix(KM.DISPLACEMENT_X)
                node.Fix(KM.DISPLACEMENT_Y)
            node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, [0., -10., 0.])
            node.SetSolutionStepValue(KM.TEMPERATURE, 30.)
            node.SetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE, 20.)
        properties = part.CreateNewProperties(1)
        properties[KM.YOUNG_MODULUS] = 10000.
        properties[KM.POISSON_RATIO] = 0.2
        properties[KM.DENSITY] = 2.
        properties[KM.THICKNESS] = 1.
        properties[Dam.THERMAL_EXPANSION] = 1.e-5
        properties[KM.CONSTITUTIVE_LAW] = (Dam.ThermalLinearElastic2DPlaneStrain()
                                           if thermal else SMA.LinearElasticPlaneStrain2DLaw())
        connectivity = ([[1, 2, 4, 3], [3, 4, 6, 5]] if quadrilaterals else
                        [[1, 2, 4], [1, 4, 3], [3, 4, 6], [3, 6, 5]])
        name = "SmallDisplacementElement2D" + ("4N" if quadrilaterals else "3N")
        for element_id, nodes in enumerate(connectivity, 1):
            part.CreateNewElement(name, element_id, nodes, properties)
        builder = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        strategy = KM.ResidualBasedLinearStrategy(part, scheme, builder, True, False, False, False)
        strategy.SetEchoLevel(0)
        strategy.Initialize()
        strategy.Check()
        return model, part, strategy

    @staticmethod
    def _values(part, variable):
        return [[list(value) for value in element.CalculateOnIntegrationPoints(variable, part.ProcessInfo)]
                for element in part.Elements]

    @staticmethod
    def _displacements(part):
        return [value for node in part.Nodes
                for value in list(node.GetSolutionStepValue(KM.DISPLACEMENT))[:2]]

    @staticmethod
    def _residual(part):
        result = [0.] * (2 * part.NumberOfNodes())
        for element in part.Elements:
            lhs, rhs = KM.Matrix(), KM.Vector()
            element.CalculateLocalSystem(lhs, rhs, part.ProcessInfo)
            for i, node in enumerate(element.GetGeometry()):
                for j in range(2):
                    result[2 * (node.Id - 1) + j] += rhs[2 * i + j]
        return result

    @staticmethod
    def _assign_states(part, stresses, strains):
        # The standard process accepts one state per model part. Give each
        # element its own submodelpart to demonstrate element-wise assignment.
        for i, element in enumerate(part.Elements):
            subpart = part.CreateSubModelPart("initial_state_{}".format(i))
            subpart.AddElement(element)
            process = KM.SetInitialStateProcess2D(
                subpart, KM.Vector(strains[i]), KM.Vector(stresses[i][0]), KM.Matrix(2, 2))
            process.ExecuteInitializeSolutionStep()

    def _assert_close(self, actual, expected, tolerance=1.e-10):
        self.assertEqual(len(actual), len(expected))
        self.assertLessEqual(max((abs(a - b) for a, b in zip(actual, expected)), default=0.), tolerance)

    @staticmethod
    def _flatten(values):
        return [component for element in values for point in element for component in point]

    def _continuation(self, thermal=False, quadrilaterals=False):
        model, source, reference = self._body(thermal=thermal, quadrilaterals=quadrilaterals)
        reference.Solve()
        stress = self._values(source, KM.CAUCHY_STRESS_VECTOR)
        strain = self._values(source, KM.GREEN_LAGRANGE_STRAIN_VECTOR)
        displacement = self._displacements(source)
        baseline_rhs = self._residual(source)
        self.assertLess(max(abs(v) for v in baseline_rhs[4:]), 1.e-10)
        self.assertAlmostEqual(baseline_rhs[1] + baseline_rhs[3], -40., delta=1.e-10)
        reactions = [value for node in source.Nodes for value in
                     list(node.GetSolutionStepValue(KM.REACTION))[:2]]
        self._assert_close(reactions[:4], [-value for value in baseline_rhs[:4]])
        # With body acceleration zero, the element RHS is minus the internal force.
        for node in source.Nodes:
            node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, [0., 0., 0.])
        minus_internal = self._residual(source)
        external = [a - b for a, b in zip(baseline_rhs, minus_internal)]
        self.assertAlmostEqual(sum(external[1::2]), -40., delta=1.e-10)
        self.assertAlmostEqual(sum(external[0::2]), 0., delta=1.e-10)
        for node in source.Nodes:
            node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, [2., -10., 0.])
        reference.Solve()
        reference_increment = [a - b for a, b in zip(self._displacements(source), displacement)]
        reference_stress = self._flatten(self._values(source, KM.CAUCHY_STRESS_VECTOR))
        reference_reactions = self._residual(source)[:4]
        for i, node in enumerate(source.Nodes):
            node.SetSolutionStepValue(KM.DISPLACEMENT, [displacement[2*i], displacement[2*i+1], 0.])
        thermal_strain = 1.e-5 * 1.2 * 10. if thermal else 0.
        modes = ("displacement", "combined") if quadrilaterals else ("displacement", "stress", "combined")
        for mode in modes:
            with self.subTest(thermal=thermal, mode=mode):
                target_model, target, strategy = self._body(thermal=thermal, quadrilaterals=quadrilaterals)
                for node in target.Nodes:
                    node.SetSolutionStepValue(KM.TEMPERATURE, 0.)
                    node.SetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE, 0.)
                mapper = Dam.MappingVariables2DUtilities()
                mapper.MappingThermalModelParts(source, target, True, True)
                for node in target.Nodes:
                    self.assertAlmostEqual(node.GetSolutionStepValue(KM.TEMPERATURE), 30.)
                    self.assertAlmostEqual(node.GetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE), 20.)
                if mode != "stress":
                    mapper.MappingMechanicalModelParts(source, target, True, False)
                if mode != "displacement":
                    offsets = []
                    for element_strain in strain:
                        # In this law the reported strain already has thermal strain subtracted.
                        base = element_strain[0] if mode == "combined" else [-thermal_strain, -thermal_strain, 0.]
                        offsets.append(base)
                    self._assign_states(target, stress, offsets)
                initial_u = self._displacements(target)
                expected_u = [0.] * len(displacement) if mode == "stress" else displacement
                self._assert_close(initial_u, expected_u)
                self._assert_close(self._flatten(self._values(target, KM.CAUCHY_STRESS_VECTOR)), self._flatten(stress))
                self._assert_close(self._residual(target), baseline_rhs)
                strategy.Solve()
                self._assert_close(self._displacements(target), initial_u)
                self._assert_close(self._flatten(self._values(target, KM.CAUCHY_STRESS_VECTOR)), self._flatten(stress))
                self._assert_close(self._residual(target), baseline_rhs)
                reactions = [value for node in target.Nodes for value in
                             list(node.GetSolutionStepValue(KM.REACTION))[:2]]
                self._assert_close(reactions[:4], [-value for value in baseline_rhs[:4]])
                for node in target.Nodes:
                    node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, [2., -10., 0.])
                strategy.Solve()
                increments = [a - b for a, b in zip(self._displacements(target), initial_u)]
                self._assert_close(increments, reference_increment)
                self._assert_close(self._flatten(self._values(target, KM.CAUCHY_STRESS_VECTOR)), reference_stress)
                self._assert_close(self._residual(target)[:4], reference_reactions)

    def test_triangle_elastic_continuation(self):
        self._continuation()

    def test_uncompensated_stress_and_displacement_double_count(self):
        model, source, strategy = self._body()
        strategy.Solve()
        stress = self._values(source, KM.CAUCHY_STRESS_VECTOR)
        target_model, target, target_strategy = self._body()
        Dam.MappingVariables2DUtilities().MappingMechanicalModelParts(source, target, True, False)
        self._assign_states(target, stress, [[0., 0., 0.] for _ in stress])
        actual = self._flatten(self._values(target, KM.CAUCHY_STRESS_VECTOR))
        self._assert_close(actual, [2. * value for value in self._flatten(stress)])
        self.assertGreater(max(abs(value) for value in self._residual(target)[4:]), 1.)

    def test_triangle_thermoelastic_continuation(self):
        self._continuation(thermal=True)

    def test_quadrilateral_displacement_and_combined_continuation(self):
        self._continuation(quadrilaterals=True)

    def test_quadrilateral_elementwise_stress_is_not_gauss_point_transfer(self):
        model, source, strategy = self._body(quadrilaterals=True)
        strategy.Solve()
        stress = self._values(source, KM.CAUCHY_STRESS_VECTOR)
        target_model, target, target_strategy = self._body(quadrilaterals=True)
        self._assign_states(target, stress, [[0., 0., 0.] for _ in stress])
        initial_stress = self._flatten(self._values(target, KM.CAUCHY_STRESS_VECTOR))
        self.assertGreater(max(abs(a - b) for a, b in zip(initial_stress, self._flatten(stress))), 1.)
        self.assertGreater(max(abs(v) for v in self._residual(target)[4:]), 1.)
        target_strategy.Solve()
        self.assertGreater(max(abs(v) for v in self._displacements(target)), 1.e-4)
        self.assertLess(max(abs(v) for v in self._residual(target)[4:]), 1.e-10)


if __name__ == "__main__":
    KratosUnittest.main()
