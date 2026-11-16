import math

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Application variables that are registered globally but are not exposed as
# module attributes.
_INITIAL_JOINT_WIDTH = KM.KratosGlobals.GetVariable("INITIAL_JOINT_WIDTH")
_NORMAL_STIFFNESS = KM.KratosGlobals.GetVariable("NORMAL_STIFFNESS")
_SHEAR_STIFFNESS = KM.KratosGlobals.GetVariable("SHEAR_STIFFNESS")
_PENALTY_STIFFNESS = KM.KratosGlobals.GetVariable("PENALTY_STIFFNESS")
_STATE_VARIABLE = KM.KratosGlobals.GetVariable("STATE_VARIABLE")
_DAMAGE_THRESHOLD = KM.KratosGlobals.GetVariable("DAMAGE_THRESHOLD")
_CRITICAL_DISPLACEMENT = KM.KratosGlobals.GetVariable("CRITICAL_DISPLACEMENT")
_YIELD_STRESS = KM.KratosGlobals.GetVariable("YIELD_STRESS")
_FRICTION_COEFFICIENT = KM.KratosGlobals.GetVariable("FRICTION_COEFFICIENT")


class TestDamUpliftJointState(KratosUnittest.TestCase):
    """Uplift pressure must follow the current mechanical aperture.

    The uplift process used to read STATE_VARIABLE from the joint elements and
    treat a value of exactly 1.0 as "broken". The Dam interface element did not
    provide that output at all, so every joint was seen as closed. On top of
    that, STATE_VARIABLE is a constitutive loading-history (damage) measure in
    the laws that expose it, not the current aperture. These tests pin the
    mechanical-aperture criterion: a joint point is open when its current
    aperture exceeds the initial joint width.
    """

    _length = 1.0
    _n_elem = 4
    _initial_width = 1.0e-3
    _opening = 1.0e-3
    _gamma = 10000.0
    _water_level = 10.0
    _full_pressure = _gamma * _water_level

    # ------------------------------------------------------------------
    # Model construction helpers
    # ------------------------------------------------------------------

    def _BuildMainModelPart(self, model):
        main = model.CreateModelPart("Uplift")
        main.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_PRESSURE)
        main.ProcessInfo[KratosDam.TIME_UNIT_CONVERTER] = 1.0

        for i in range(self._n_elem + 1):
            x = self._length * i / self._n_elem
            main.CreateNewNode(100 + i, x, 0.0, 0.0)

        return main

    def _CreateJointProperties(self, model_part, constitutive_law):
        props = model_part.CreateNewProperties(1)
        props.SetValue(KM.CONSTITUTIVE_LAW, constitutive_law)
        props.SetValue(_INITIAL_JOINT_WIDTH, self._initial_width)
        props.SetValue(_NORMAL_STIFFNESS, 3.0e5)
        props.SetValue(_SHEAR_STIFFNESS, 3.0e5)
        props.SetValue(_PENALTY_STIFFNESS, 1.0e3)
        props.SetValue(KM.DENSITY, 0.0)
        props.SetValue(KM.THICKNESS, 1.0)
        return props

    def _BuildJoint2D(self, model, constitutive_law):
        joint = model.CreateModelPart("Joint")
        joint.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        lower = []
        upper = []
        for i in range(self._n_elem + 1):
            x = self._length * i / self._n_elem
            node_lower = joint.CreateNewNode(2 * i + 1, x, 0.0, 0.0)
            node_upper = joint.CreateNewNode(2 * i + 2, x, 0.0, 0.0)
            for node in (node_lower, node_upper):
                node.AddDof(KM.DISPLACEMENT_X)
                node.AddDof(KM.DISPLACEMENT_Y)
            lower.append(node_lower)
            upper.append(node_upper)

        props = self._CreateJointProperties(joint, constitutive_law)
        for e in range(self._n_elem):
            joint.CreateNewElement(
                "SmallDisplacementInterfaceElement2D4N",
                e + 1,
                [
                    lower[e].Id,
                    lower[e + 1].Id,
                    upper[e + 1].Id,
                    upper[e].Id,
                ],
                props,
            )

        for element in joint.Elements:
            element.Initialize(joint.ProcessInfo)

        return joint, lower, upper

    def _BuildJoint3D(self, model, constitutive_law):
        joint = model.CreateModelPart("Joint3D")
        joint.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        # A single 3D8N interface element on the plane z=0. Nodes 0..3 are the
        # lower face and nodes 4..7 the opposite face; the local normal is +z.
        coordinates = [
            (0.0, 0.0, 0.0),
            (self._length, 0.0, 0.0),
            (self._length, 1.0, 0.0),
            (0.0, 1.0, 0.0),
        ]
        nodes = []
        for i, (x, y, z) in enumerate(coordinates + coordinates):
            node = joint.CreateNewNode(i + 1, x, y, z)
            for dof in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
                node.AddDof(dof)
            nodes.append(node)

        props = self._CreateJointProperties(joint, constitutive_law)
        props.SetValue(_DAMAGE_THRESHOLD, 0.2)
        props.SetValue(_CRITICAL_DISPLACEMENT, 0.1)
        props.SetValue(KM.YOUNG_MODULUS, 1.0e6)
        props.SetValue(_YIELD_STRESS, 1.0e5)
        props.SetValue(_FRICTION_COEFFICIENT, 0.0)

        joint.CreateNewElement(
            "SmallDisplacementInterfaceElement3D8N",
            1,
            [node.Id for node in nodes],
            props,
        )
        for element in joint.Elements:
            element.Initialize(joint.ProcessInfo)

        return joint, nodes

    # ------------------------------------------------------------------
    # Prescribed kinematics and reference solution
    # ------------------------------------------------------------------

    @staticmethod
    def _SetVerticalOpening(upper, x_front, opening):
        for node in upper:
            displacement = KM.Array3(0.0)
            if node.X <= x_front + 1.0e-12:
                displacement[1] = opening
            node.SetSolutionStepValue(KM.DISPLACEMENT, displacement)

    @staticmethod
    def _SetOutOfPlaneOpening(nodes, opening):
        for index, node in enumerate(nodes):
            displacement = KM.Array3(0.0)
            if index >= 4:
                displacement[2] = opening
            node.SetSolutionStepValue(KM.DISPLACEMENT, displacement)

    def _CreateProcess(self, main, joint):
        settings = KM.Parameters(
            """
            {
                "model_part_name"                                  : "Uplift",
                "variable_name"                                    : "POSITIVE_FACE_PRESSURE",
                "joint_group_name"                                 : "Joint",
                "Gravity_Direction"                                : "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction" : 0.0,
                "Upstream_Coordinate"                              : [0.0, 0.0, 0.0],
                "Downstream_Coordinate"                            : [%.15g, 0.0, 0.0],
                "Upstream_Longitudinal_Coordinate"                 : [0.0, 1.0, 0.0],
                "Spe_weight"                                       : %.15g,
                "Water_level"                                      : %.15g,
                "Water_Table"                                      : 0,
                "Drains"                                           : false
            }
            """
            % (self._length, self._gamma, self._water_level)
        )
        return KratosDam.DamUpliftConditionLoadProcess(main, joint, settings)

    def _ExpectedPressure(self, x, x_front):
        """Hydrostatic pressure over the open length and linear decay over the
        remaining ligament; full hydrostatic when the joint is fully open."""
        if x_front >= self._length:
            return self._full_pressure
        if x < x_front:
            return self._full_pressure
        return self._full_pressure * (
            1.0 - (x - x_front) / (self._length - x_front)
        )

    @staticmethod
    def _TrapezoidalResultant(model_part):
        xs = []
        pressures = []
        for node in model_part.Nodes:
            xs.append(node.X)
            pressures.append(
                node.GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE)
            )
        return sum(
            0.5 * (xs[i + 1] - xs[i]) * (pressures[i + 1] + pressures[i])
            for i in range(len(xs) - 1)
        )

    def _CheckAnalyticalDistribution(self, main, x_front):
        for node in main.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE),
                self._ExpectedPressure(node.X, x_front),
                delta=1.0e-9,
            )

        expected_resultant = self._full_pressure * (
            self._length + x_front
        ) / 2.0
        if x_front >= self._length:
            expected_resultant = self._full_pressure * self._length
        self.assertAlmostEqual(
            self._TrapezoidalResultant(main), expected_resultant, delta=1.0e-6
        )

    # ------------------------------------------------------------------
    # Tests
    # ------------------------------------------------------------------

    def test_closed_joint_keeps_linear_uplift(self):
        """Test A: a mechanically closed joint keeps the uncracked pressure."""
        model = KM.Model()
        main = self._BuildMainModelPart(model)
        joint, _, upper = self._BuildJoint2D(model, KratosPoro.ElasticCohesive2DLaw())

        self._SetVerticalOpening(upper, -1.0, self._opening)

        process = self._CreateProcess(main, joint)
        process.ExecuteBeforeSolutionLoop()

        # Closed junction: full pressure at the upstream face decaying linearly
        # to zero at the downstream toe. It is not zero.
        self._CheckAnalyticalDistribution(main, 0.0)
        for node in main.Nodes:
            self.assertFalse(
                math.isnan(node.GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE))
            )

    def test_partially_open_joint_uses_analytical_crack_front(self):
        """Test B: opened integration points define the crack front."""
        model = KM.Model()
        main = self._BuildMainModelPart(model)
        joint, _, upper = self._BuildJoint2D(model, KratosPoro.ElasticCohesive2DLaw())

        crack_length = 0.5
        self._SetVerticalOpening(upper, crack_length, self._opening)

        process = self._CreateProcess(main, joint)
        process.ExecuteBeforeSolutionLoop()

        self._CheckAnalyticalDistribution(main, crack_length)

        # The process must have detected the opened half, not a closed joint.
        self.assertAlmostEqual(
            main.GetNode(100).GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE),
            self._full_pressure,
            delta=1.0e-9,
        )
        self.assertAlmostEqual(
            main.GetNode(103).GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE),
            self._full_pressure * 0.5,
            delta=1.0e-9,
        )

    def test_fully_open_joint_gives_full_hydrostatic_pressure(self):
        """Test C: a joint opened up to the toe carries full pressure."""
        model = KM.Model()
        main = self._BuildMainModelPart(model)
        joint, _, upper = self._BuildJoint2D(model, KratosPoro.ElasticCohesive2DLaw())

        self._SetVerticalOpening(upper, self._length, self._opening)

        process = self._CreateProcess(main, joint)
        process.ExecuteBeforeSolutionLoop()

        for node in main.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE),
                self._full_pressure,
                delta=1.0e-9,
            )
        self.assertAlmostEqual(
            self._TrapezoidalResultant(main),
            self._full_pressure * self._length,
            delta=1.0e-6,
        )

    def test_reopening_and_closure_follow_the_current_aperture(self):
        """Test D1: the state is recomputed at every solution step."""
        model = KM.Model()
        main = self._BuildMainModelPart(model)
        joint, _, upper = self._BuildJoint2D(model, KratosPoro.ElasticCohesive2DLaw())

        process = self._CreateProcess(main, joint)

        self._SetVerticalOpening(upper, 0.5, self._opening)
        process.ExecuteBeforeSolutionLoop()
        self._CheckAnalyticalDistribution(main, 0.5)

        self._SetVerticalOpening(upper, -1.0, self._opening)
        process.ExecuteInitializeSolutionStep()
        self._CheckAnalyticalDistribution(main, 0.0)

        self._SetVerticalOpening(upper, 0.5, self._opening)
        process.ExecuteInitializeSolutionStep()
        self._CheckAnalyticalDistribution(main, 0.5)

    def test_damage_is_not_used_as_opening_indicator(self):
        """Test D2: a damaged but mechanically closed joint is not broken.

        STATE_VARIABLE is a monotonic constitutive loading-history measure in
        the laws that provide it. The uplift process must follow the current
        aperture and not that history.
        """
        model = KM.Model()
        main = self._BuildMainModelPart(model)
        joint, nodes = self._BuildJoint3D(model, KratosPoro.BilinearCohesive3DLaw())
        element = joint.GetElement(1)

        # The binding validates the values against the geometry output points,
        # while the element stores one value per constitutive-law point.
        n_output_points = len(element.GetIntegrationPoints())
        element.SetValuesOnIntegrationPoints(
            _STATE_VARIABLE, [1.0] * n_output_points, joint.ProcessInfo
        )
        states = list(
            element.CalculateOnIntegrationPoints(_STATE_VARIABLE, joint.ProcessInfo)
        )
        self.assertTrue(states)
        for state in states:
            self.assertAlmostEqual(state, 1.0, delta=1.0e-15)

        process = self._CreateProcess(main, joint)

        # Fully damaged, but mechanically closed.
        self._SetOutOfPlaneOpening(nodes, 0.0)
        process.ExecuteBeforeSolutionLoop()
        self._CheckAnalyticalDistribution(main, 0.0)

        # Same damage state, now mechanically open.
        self._SetOutOfPlaneOpening(nodes, self._opening)
        process.ExecuteInitializeSolutionStep()
        for node in main.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE),
                self._full_pressure,
                delta=1.0e-9,
            )

        # The constitutive history has not been modified by the process.
        states = list(
            element.CalculateOnIntegrationPoints(_STATE_VARIABLE, joint.ProcessInfo)
        )
        self.assertTrue(states)
        for state in states:
            self.assertAlmostEqual(state, 1.0, delta=1.0e-15)

    def test_uplift_detection_is_constitutive_law_independent(self):
        """Test E: laws without STATE_VARIABLE give the same uplift."""
        pressure_sets = []
        for law in (
            KratosPoro.ElasticCohesive2DLaw(),
            KratosPoro.IsotropicDamageCohesive2DLaw(),
        ):
            model = KM.Model()
            main = self._BuildMainModelPart(model)
            joint, _, upper = self._BuildJoint2D(model, law)
            self._SetVerticalOpening(upper, 0.5, self._opening)

            process = self._CreateProcess(main, joint)
            process.ExecuteBeforeSolutionLoop()

            self._CheckAnalyticalDistribution(main, 0.5)
            pressure_sets.append(
                [
                    node.GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE)
                    for node in main.Nodes
                ]
            )

        for first, second in zip(pressure_sets[0], pressure_sets[1]):
            self.assertAlmostEqual(first, second, delta=1.0e-12)

    def test_element_reports_mechanical_aperture_per_integration_point(self):
        """Test F: JOINT_WIDTH is the current mechanical aperture."""
        model = KM.Model()
        main = self._BuildMainModelPart(model)
        joint, _, upper = self._BuildJoint2D(model, KratosPoro.ElasticCohesive2DLaw())

        first_element = joint.GetElement(1)

        self._SetVerticalOpening(upper, -1.0, self._opening)
        widths = first_element.CalculateOnIntegrationPoints(
            KratosPoro.JOINT_WIDTH, joint.ProcessInfo
        )
        self.assertEqual(len(widths), 2)
        for width in widths:
            self.assertAlmostEqual(width, self._initial_width, delta=1.0e-15)

        self._SetVerticalOpening(upper, self._length, self._opening)
        widths = first_element.CalculateOnIntegrationPoints(
            KratosPoro.JOINT_WIDTH, joint.ProcessInfo
        )
        for width in widths:
            self.assertAlmostEqual(
                width, self._initial_width + self._opening, delta=1.0e-15
            )

        # The elastic cohesive law does not provide STATE_VARIABLE; the element
        # reports zero instead of leaking an uninitialized value.
        states = first_element.CalculateOnIntegrationPoints(
            _STATE_VARIABLE, joint.ProcessInfo
        )
        self.assertEqual(list(states), [0.0, 0.0])


if __name__ == "__main__":
    KratosUnittest.main()
