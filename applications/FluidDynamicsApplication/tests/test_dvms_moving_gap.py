import copy

import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as CFD
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication import dvms_gap_history
from KratosMultiphysics.FluidDynamicsApplication.two_fluid_interface_time_coupling import (
    PredictorContext, TwoFluidInterfaceTimeCoupling)
from KratosMultiphysics.FluidDynamicsApplication.two_fluid_level_set_geometry import extract_height_contour
from test_one_phase_moving_gap import assemble_residual, make_lifecycle


class TestDVMSMovingGap(KratosUnittest.TestCase):
    def test_one_and_two_phase_geometry(self):
        for phase in ("one_phase", "two_phase"):
            lifecycle = make_lifecycle(phase, fluid_element="dvms")
            self.assertTrue(all(CFD.IsTimeIntegratedDVMS2D4N(e) for e in lifecycle.computing.Elements))
            gap = lifecycle.root.GetSubModelPart("TwoFluidGapCoupling")
            self.assertGreater(gap.GetSubModelPart("GapVolumes").NumberOfElements(), 0)
            self.assertGreater(gap.GetSubModelPart("GapInterfaces").NumberOfConditions(), 0)
            self.assertEqual(gap.GetSubModelPart("PhysicalInterfaces").NumberOfConditions() == 0,
                             phase == "one_phase")
        self.assertEqual(make_lifecycle().fluid_element, "qsvms")

    def test_hydrostatic_and_zero_feedback(self):
        lifecycle = make_lifecycle(mass_source=True, fluid_element="dvms", mass_source_strategy="relaxed")
        lifecycle._prepare_mass_source()
        self.assertEqual(lifecycle._pending_mass_source["mass_source_water_rate"], 0.0)
        self.assertLess(assemble_residual(lifecycle.computing), 1e-8)

    def test_transfer_preserves_constant_and_phase_identity(self):
        lifecycle = make_lifecycle("two_phase", fluid_element="dvms")
        water_ids = {node.Id for node in lifecycle.water.Nodes}
        for element in lifecycle.computing.Elements:
            water = next(iter(element.GetGeometry())).Id in water_ids
            state = CFD.GetDVMSAcceptedSubscaleState(element)
            values = KM.Matrix(state.Size1(), 3)
            for i in range(state.Size1()):
                values[i, 0] = 2.0 if water else -3.0
                values[i, 1] = -1.0 if water else 4.0
                values[i, 2] = 0.25 if water else -0.5
            CFD.SetDVMSAcceptedSubscaleState(element, values, 7)
        snapshot = dvms_gap_history.capture(lifecycle.computing, lifecycle.water, lifecycle.air)
        old_ids = {node.Id for node in lifecycle.phases[0].Nodes}
        lifecycle._build_geometry([(0.0, 0.38), (1.0, 0.38)])
        dvms_gap_history.restore(lifecycle.computing, lifecycle.water, lifecycle.air, snapshot)
        self.assertNotEqual(old_ids, {node.Id for node in lifecycle.phases[0].Nodes})
        for element in lifecycle.computing.Elements:
            water = next(iter(element.GetGeometry())).Id in water_ids
            state = CFD.GetDVMSAcceptedSubscaleState(element)
            self.assertEqual(CFD.GetDVMSLastCommittedSubscaleStep(element), 7)
            for i in range(state.Size1()):
                self.assertAlmostEqual(state[i, 3], 2.0 if water else -3.0)
                self.assertAlmostEqual(state[i, 4], -1.0 if water else 4.0)
                self.assertAlmostEqual(state[i, 5], 0.25 if water else -0.5)

    def test_affine_fit_and_singular_fallback(self):
        samples = [[x, y, 1.0, 1.0 + 2*x - y, -2.0 + x + 3*y, 0.0]
                   for x, y in ((0., 0.), (1., 0.), (0., 1.), (1., 1.))]
        velocity = dvms_gap_history._solve_affine(samples, 0.3, 0.2)
        self.assertAlmostEqual(velocity[0], 1.4)
        self.assertAlmostEqual(velocity[1], -1.1)
        self.assertIsNone(dvms_gap_history._solve_affine(samples[:2], 0.3, 0.2))

    def test_checkpoint_and_isolation(self):
        accepted = make_lifecycle(fluid_element="dvms")
        predictor = make_lifecycle(fluid_element="dvms")
        state = accepted.export_lifecycle_state()
        predictor.restore_lifecycle_state(copy.deepcopy(state))
        predictor.prepare_fixed_geometry_predictor(0.1)
        self.assertEqual(accepted.export_lifecycle_state(), state)
        bad = copy.deepcopy(state)
        bad["dvms_subscale_history"] = None
        with self.assertRaisesRegex(ValueError, "missing accepted subscale"):
            make_lifecycle(fluid_element="dvms").restore_lifecycle_state(bad)
        with self.assertRaisesRegex(ValueError, "checkpoint"):
            make_lifecycle().restore_lifecycle_state(state)

    def test_predicted_velocity_full_step_is_independent(self):
        self._check_predicted_velocity_full_step("dvms")
        self._check_predicted_velocity_full_step("qsvms")

    def _check_predicted_velocity_full_step(self, fluid_element):
        accepted = make_lifecycle(mass_source=True, fluid_element=fluid_element, mass_source_strategy="relaxed")
        predictor = make_lifecycle(mass_source=True, fluid_element=fluid_element, mass_source_strategy="relaxed")
        for lifecycle in (accepted, predictor):
            for node in lifecycle.water.Nodes:
                for slot in range(3):
                    node.SetSolutionStepValue(KM.BODY_FORCE, slot, [0.1, -9.81, 0.0])

        def advance(time, sampler):
            part = accepted.transport.part
            changes = {node.Id: sampler(node.X, node.Y, int(node.GetSolutionStepValue(KM.DISTANCE) > 0))[1]
                       for node in part.Nodes}
            part.CloneTimeStep(time)
            for node in part.Nodes:
                node.SetSolutionStepValue(KM.DISTANCE,
                                          node.GetSolutionStepValue(KM.DISTANCE) - 0.1 * changes[node.Id])
            return extract_height_contour(part, accepted.tolerance)

        accepted.transport.advance = advance
        settings = KM.Parameters('''{"relative_tolerance":1e-8,"absolute_tolerance":1e-10,
            "maximum_iterations":50,"maximum_free_dof_residual":1e-7,"echo_level":0}''')
        before = accepted.export_lifecycle_state()

        def factory(instance):
            self.assertEqual(instance.export_lifecycle_state(), before)
            return PredictorContext(predictor, KM.SkylineLUFactorizationSolver(), settings.Clone())

        coupling = TwoFluidInterfaceTimeCoupling(accepted, "predicted_velocity", factory)
        residual, _ = coupling.advance(0.1, KM.SkylineLUFactorizationSolver(), settings)
        self.assertLess(residual, 1e-7)
        if fluid_element == "dvms":
            self.assertEqual(accepted._dvms_snapshot["step"], 1)
            self.assertEqual(predictor._dvms_snapshot["step"], 1)
            self.assertIsNot(accepted._dvms_snapshot, predictor._dvms_snapshot)
        self.assertEqual(coupling.diagnostics["fluid_solves"], 2)

    def test_continuity_source_integral(self):
        lifecycle = make_lifecycle(fluid_element="dvms")
        element = next(iter(lifecycle.root.GetSubModelPart("TwoFluidGapCoupling")
                            .GetSubModelPart("GapVolumes").Elements))
        lhs, before, after = KM.Matrix(), KM.Vector(), KM.Vector()
        element.CalculateLocalSystem(lhs, before, lifecycle.root.ProcessInfo)
        element.SetValue(CFD.CONTINUITY_SOURCE, 2.0)
        element.CalculateLocalSystem(lhs, after, lifecycle.root.ProcessInfo)
        state = CFD.GetDVMSAcceptedSubscaleState(element)
        area = sum(state[i, 2] for i in range(state.Size1()))
        self.assertAlmostEqual(sum(after[3*i+2] - before[3*i+2] for i in range(4)), 2.0*area)
        self.assertGreater(sum(abs(after[3*i+d] - before[3*i+d]) for i in range(4) for d in range(2)), 0.0)
        element.SetValue(CFD.CONTINUITY_SOURCE, 0.0)
        previous = KM.Matrix(state.Size1(), 3)
        for i in range(state.Size1()):
            previous[i, 0], previous[i, 1], previous[i, 2] = state[i, 3], state[i, 4], 2.0
        CFD.SetDVMSAcceptedSubscaleState(element, previous, -1)
        element.CalculateLocalSystem(lhs, after, lifecycle.root.ProcessInfo)
        self.assertAlmostEqual(sum(abs(after[3*i+2] - before[3*i+2]) for i in range(4)), 0.0)
        self.assertGreater(sum(abs(after[3*i+d] - before[3*i+d]) for i in range(4) for d in range(2)), 0.0)

    def test_source_power_preserves_dvms_state(self):
        lifecycle = make_lifecycle("two_phase", mass_source=True, fluid_element="dvms",
                                   mass_source_strategy="relaxed")
        gap = lifecycle.root.GetSubModelPart("TwoFluidGapCoupling").GetSubModelPart("GapVolumes")
        for element in gap.Elements:
            element.SetValue(CFD.CONTINUITY_SOURCE, 0.25 if element.Properties.Id == 1 else -0.25)
        before = dvms_gap_history.capture(lifecycle.computing, lifecycle.water, lifecycle.air)
        expected = 0.0
        for element in lifecycle.computing.Elements:
            rhs = CFD.CalculateDVMSContinuitySourceRHS(element, lifecycle.computing.ProcessInfo)
            for i, node in enumerate(element.GetGeometry()):
                velocity = node.GetSolutionStepValue(KM.VELOCITY)
                expected += sum(velocity[d] * rhs[3*i+d] for d in range(2))
                expected += node.GetSolutionStepValue(KM.PRESSURE) * rhs[3*i+2]
        self.assertAlmostEqual(lifecycle._source_power(), expected)
        self.assertEqual(before, dvms_gap_history.capture(lifecycle.computing, lifecycle.water, lifecycle.air))


if __name__ == "__main__":
    KratosUnittest.main()
