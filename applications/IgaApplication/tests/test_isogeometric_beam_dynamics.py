"""Analytical vibration benchmarks for the simplified IGA beam inertia."""
import math

import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestIsogeometricBeamDynamics(KratosUnittest.TestCase):
    def run_constrained_vibration(self, steps_per_period, motion="axial"):
        model = KM.Model()
        beam = model.CreateModelPart("beam", 2)
        for variable in (KM.DISPLACEMENT, KM.ROTATION, KM.VELOCITY, KM.ACCELERATION,
                         KM.ANGULAR_VELOCITY, KM.ANGULAR_ACCELERATION):
            beam.AddNodalSolutionStepVariable(variable)
        beam.ProcessInfo[KM.DOMAIN_SIZE] = 3
        length, young, density, area = 2.0, 1000.0, 3.0, 2.0
        amplitude = 1e-7  # Linear limit of the geometrically nonlinear element.
        points = KM.NodesVector()
        for i in range(3):
            node = beam.CreateNewNode(i + 1, length * i / 2, 0.0, 0.0)
            points.append(node)
            for variable in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z, KM.ROTATION_X):
                node.AddDof(variable)
                node.Fix(variable)
        tip = beam.GetNode(3)
        dof_variable = {"axial": KM.DISPLACEMENT_X, "bending": KM.DISPLACEMENT_Y,
                    "torsion": KM.ROTATION_X}[motion]
        velocity_variable = KM.ANGULAR_VELOCITY_X if motion == "torsion" else (
            KM.VELOCITY_Y if motion == "bending" else KM.VELOCITY_X)
        acceleration_variable = KM.ANGULAR_ACCELERATION_X if motion == "torsion" else (
            KM.ACCELERATION_Y if motion == "bending" else KM.ACCELERATION_X)
        tip.Free(dof_variable)
        curve = KM.NurbsCurveGeometry3D(points, 2, KM.Vector([0.0, 0.0, 1.0, 1.0]))
        properties = beam.CreateNewProperties(1)
        for variable, value in ((KM.YOUNG_MODULUS, young), (KM.POISSON_RATIO, 0.25),
                                (KM.DENSITY, density), (IGA.CROSS_AREA, area),
                                (IGA.I_N, 0.4), (IGA.I_V, 0.7), (IGA.I_T, 0.5)):
            properties.SetValue(variable, value)
        properties.SetValue(IGA.BEAM_MASS_FORMULATION, "simplified")
        properties.SetValue(KM.CONSTITUTIVE_LAW, IGA.BernoulliBeamElasticConstitutiveLaw())
        properties.SetValue(IGA.T_0, KM.Vector([1.0, 0.0, 0.0]))
        properties.SetValue(IGA.N_0, KM.Vector([0.0, 0.0, 1.0]))
        orientation = KM.Matrix(2, 4, 0.0)
        orientation[1, 0] = 1.0
        orientation[0, 3] = orientation[1, 3] = 1.0
        properties.SetValue(IGA.LOCAL_AXIS_ORIENTATION, orientation)
        quadrature = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(quadrature, 3)
        for i in range(len(quadrature)):
            beam.CreateNewElement("IsogeometricBeamElement", i + 1, quadrature[i], properties)
        # Only the selected tip coefficient is free: N_tip=(s/L)^2.
        mass = density * area * length / 5
        stiffness = 4 * young * area / (3 * length)
        if motion == "bending":
            stiffness = 4 * young * 0.4 / length**3
        elif motion == "torsion":
            mass = density * (0.4 + 0.7) * length / 5
            shear = young / (2 * (1 + 0.25))
            stiffness = 4 * shear * 0.5 / (3 * length)
        omega = math.sqrt(stiffness / mass)
        period = 2 * math.pi / omega
        dt = period / steps_per_period
        for step in (0, 1):
            tip.SetSolutionStepValue(dof_variable, step, amplitude)
            tip.SetSolutionStepValue(acceleration_variable, step, -omega**2 * amplitude)
        scheme = KM.ResidualBasedBossakDisplacementScheme(0.0)  # Average-acceleration Newmark.
        criteria = KM.ResidualCriteria(1e-10, 1e-13)
        criteria.SetEchoLevel(0)
        solver = KM.ResidualBasedNewtonRaphsonStrategy(
            beam, scheme, KM.SkylineLUFactorizationSolver(), criteria, 15, False, False, False)
        solver.SetEchoLevel(0)
        if motion == "torsion":
            for element in beam.Elements:
                element.Initialize(beam.ProcessInfo)
        maximum_error = 0.0
        maximum_energy_error = 0.0
        initial_energy = 0.5 * stiffness * amplitude**2
        for step in range(1, 2 * steps_per_period + 1):
            time = step * dt
            beam.CloneTimeStep(time)
            if motion == "torsion":
                # Test-only scalar Newmark driver: the displacement Bossak scheme
                # does not update scalar twist rates/accelerations.
                previous_q = tip.GetSolutionStepValue(dof_variable, 1)
                previous_v = tip.GetSolutionStepValue(velocity_variable, 1)
                previous_a = tip.GetSolutionStepValue(acceleration_variable, 1)
                predictor = previous_q + dt * previous_v + 0.25 * dt**2 * previous_a
                q = predictor
                for iteration in range(10):
                    acceleration = 4 * (q - predictor) / dt**2
                    velocity = previous_v + 0.5 * dt * (previous_a + acceleration)
                    tip.SetSolutionStepValue(dof_variable, q)
                    tip.SetSolutionStepValue(velocity_variable, velocity)
                    tip.SetSolutionStepValue(acceleration_variable, acceleration)
                    residual, tangent = 0.0, 0.0
                    for element in beam.Elements:
                        lhs, rhs, element_mass = KM.Matrix(), KM.Vector(), KM.Matrix()
                        element.CalculateLocalSystem(lhs, rhs, beam.ProcessInfo)
                        element.CalculateMassMatrix(element_mass, beam.ProcessInfo)
                        residual += rhs[11] - element_mass[11, 11] * acceleration
                        tangent += lhs[11, 11] + 4 * element_mass[11, 11] / dt**2
                    if abs(residual) < 1e-13:
                        break
                    q += residual / tangent
                else:
                    self.fail("Scalar torsion Newmark iteration did not converge")
            else:
                solver.Solve()
            displacement = tip.GetSolutionStepValue(dof_variable)
            velocity = tip.GetSolutionStepValue(velocity_variable)
            maximum_error = max(maximum_error, abs(displacement / amplitude - math.cos(omega * time)))
            energy = 0.5 * mass * velocity**2 + 0.5 * stiffness * displacement**2
            maximum_energy_error = max(maximum_energy_error, abs(energy / initial_energy - 1))
        solver.Clear()
        return maximum_error, maximum_energy_error

    def test_simplified_mass_bending_free_vibration(self):
        coarse, coarse_energy = self.run_constrained_vibration(40, "bending")
        fine, fine_energy = self.run_constrained_vibration(80, "bending")
        self.assertLess(fine, 0.007)
        self.assertLess(fine, 0.3 * coarse)
        self.assertLess(max(coarse_energy, fine_energy), 1e-5)

    def test_simplified_mass_torsion_free_vibration(self):
        coarse, coarse_energy = self.run_constrained_vibration(40, "torsion")
        fine, fine_energy = self.run_constrained_vibration(80, "torsion")
        self.assertLess(fine, 0.007)
        self.assertLess(fine, 0.3 * coarse)
        self.assertLess(max(coarse_energy, fine_energy), 1e-5)

    def run_distributed_axial_vibration(self, spans, steps_per_period):
        import numpy as np

        # Full open knot vector for the independent basis evaluation below.
        knots = [0.0] * 3 + [i / spans for i in range(1, spans)] + [1.0] * 3
        count = spans + 2
        greville = [(knots[i + 1] + knots[i + 2]) / 2 for i in range(count)]

        def basis(x):
            if x == 1.0:
                return np.array([0.0] * (count - 1) + [1.0])
            values = np.array([float(knots[i] <= x < knots[i + 1])
                               for i in range(len(knots) - 1)])
            for degree in (1, 2):
                result = np.zeros(len(values) - 1)
                for i in range(len(result)):
                    left = knots[i + degree] - knots[i]
                    right = knots[i + degree + 1] - knots[i + 1]
                    if left > 0:
                        result[i] += (x - knots[i]) / left * values[i]
                    if right > 0:
                        result[i] += (knots[i + degree + 1] - x) / right * values[i + 1]
                values = result
            return values

        model = KM.Model()
        beam = model.CreateModelPart("beam", 2)
        for variable in (KM.DISPLACEMENT, KM.ROTATION, KM.VELOCITY, KM.ACCELERATION,
                         KM.ANGULAR_VELOCITY, KM.ANGULAR_ACCELERATION):
            beam.AddNodalSolutionStepVariable(variable)
        beam.ProcessInfo[KM.DOMAIN_SIZE] = 3
        length, young, density, area = 2.0, 1000.0, 3.0, 2.0
        amplitude = 1e-7  # Linear limit of the geometrically nonlinear element.
        points = KM.NodesVector()
        for i, coordinate in enumerate(greville):
            node = beam.CreateNewNode(i + 1, length * coordinate, 0.0, 0.0)
            points.append(node)
            for variable in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z, KM.ROTATION_X):
                node.AddDof(variable)
                node.Fix(variable)
        for node in beam.Nodes:
            if node.Id != 1:
                node.Free(KM.DISPLACEMENT_X)
        curve = KM.NurbsCurveGeometry3D(points, 2, KM.Vector(knots[1:-1]))
        properties = beam.CreateNewProperties(1)
        for variable, value in ((KM.YOUNG_MODULUS, young), (KM.POISSON_RATIO, 0.25),
                                (KM.DENSITY, density), (IGA.CROSS_AREA, area),
                                (IGA.I_N, 0.4), (IGA.I_V, 0.7), (IGA.I_T, 0.5)):
            properties.SetValue(variable, value)
        properties.SetValue(IGA.BEAM_MASS_FORMULATION, "simplified")
        properties.SetValue(KM.CONSTITUTIVE_LAW, IGA.BernoulliBeamElasticConstitutiveLaw())
        properties.SetValue(IGA.T_0, KM.Vector([1.0, 0.0, 0.0]))
        properties.SetValue(IGA.N_0, KM.Vector([0.0, 0.0, 1.0]))
        orientation = KM.Matrix(2, 4, 0.0)
        orientation[1, 0] = 1.0
        orientation[0, 3] = orientation[1, 3] = 1.0
        properties.SetValue(IGA.LOCAL_AXIS_ORIENTATION, orientation)
        quadrature = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(quadrature, 3)
        for i in range(len(quadrature)):
            beam.CreateNewElement("IsogeometricBeamElement", i + 1, quadrature[i], properties)
        # Two exact modes of rho*u_tt=E*u_ss, u(0)=0, u_s(L)=0.
        # Interpolate their initial displacement in the quadratic spline space.
        collocation = np.array([basis(x) for x in greville])
        mode_one = np.linalg.solve(collocation, np.sin(0.5 * math.pi * np.array(greville)))
        mode_two = np.linalg.solve(collocation, 0.3 * np.sin(1.5 * math.pi * np.array(greville)))
        omega = math.pi / (2 * length) * math.sqrt(young / density)
        for i, node in enumerate(beam.Nodes):
            for step in (0, 1):
                node.SetSolutionStepValue(KM.DISPLACEMENT_X, step, amplitude * (mode_one[i] + mode_two[i]))
                node.SetSolutionStepValue(KM.ACCELERATION_X, step,
                                          -amplitude * omega**2 * (mode_one[i] + 9 * mode_two[i]))
        scheme = KM.ResidualBasedBossakDisplacementScheme(0.0)
        criteria = KM.ResidualCriteria(1e-10, 1e-13)
        criteria.SetEchoLevel(0)
        solver = KM.ResidualBasedNewtonRaphsonStrategy(
            beam, scheme, KM.SkylineLUFactorizationSolver(), criteria, 15, False, False, False)
        solver.SetEchoLevel(0)
        sample_coordinates = np.linspace(0, 1, 41)
        interpolation = np.array([basis(x) for x in sample_coordinates])
        first = np.sin(0.5 * math.pi * sample_coordinates)
        second = 0.3 * np.sin(1.5 * math.pi * sample_coordinates)
        error = 0.0
        for step in range(1, steps_per_period + 1):
            time = step * 2 * math.pi / (omega * steps_per_period)
            beam.CloneTimeStep(time)
            solver.Solve()
            coefficients = np.array([node.GetSolutionStepValue(KM.DISPLACEMENT_X) / amplitude
                                     for node in beam.Nodes])
            computed = interpolation @ coefficients
            exact = first * math.cos(omega * time) + second * math.cos(3 * omega * time)
            error = max(error, float(np.max(np.abs(computed - exact))))
        solver.Clear()
        return error

    def test_simplified_mass_distributed_axial_vibration(self):
        coarse = self.run_distributed_axial_vibration(8, 200)
        fine = self.run_distributed_axial_vibration(16, 400)
        self.assertLess(fine, 0.005)
        self.assertLess(fine, 0.4 * coarse)

    def test_simplified_mass_axial_free_vibration(self):
        coarse_error, coarse_energy = self.run_constrained_vibration(40)
        fine_error, fine_energy = self.run_constrained_vibration(80)
        self.assertLess(fine_error, 0.007)
        self.assertLess(fine_error, 0.3 * coarse_error)
        self.assertLess(coarse_energy, 1e-5)
        self.assertLess(fine_energy, 1e-5)


if __name__ == "__main__":
    KratosUnittest.main()
