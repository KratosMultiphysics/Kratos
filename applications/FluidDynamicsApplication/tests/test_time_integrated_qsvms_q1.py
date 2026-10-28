#     |  /           |
#     ' /   __| _` | __|  _ \   __|
#     . \  |   (   | |   (   |\__ `
#    _|\_\_|  \__,_|\__|\___/ ____/
#                    Multi-Physics
#
#   Fluid Dynamics Application
#
#   License:         BSD License
#                    Kratos default license: kratos/license.txt
#
#   Main author:     Nicolò Antonelli
#

"""Fixed-mesh Q1 BDF2 QSVMS manufactured Navier–Stokes verification.

The existing QSVMS matrix is a Picard linearization: convection and tau are
frozen. The nonlinear strategy iterates this matrix to residual convergence.
"""
import math
import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as CFD
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication.moving_gap_navier_stokes_solver import create_qsvms_picard_strategy

MU = 0.1


def fields(x, y, time, nonlinear=False):
    a = math.exp(time)
    if nonlinear:
        base = (2*x*x*y, -2*x*y*y)
        convection = (4*x**3*y*y, 4*x*x*y**3)
        spatial_force = (2*x-4*MU*y, 2*y+4*MU*x)
        p = x*x+y*y
    else:
        base = (1+x, -y)
        convection = (1+x, y)
        spatial_force = (2.0, -3.0)
        p = 1+2*x-3*y
    velocity = [a*base[0], a*base[1], 0.0]
    # rho=1, so BODY_FORCE is acceleration: du/dt+(u.grad)u-div(sigma).
    force = [a*(base[d]+spatial_force[d])+a*a*convection[d] for d in range(2)]+[0.0]
    return velocity, a*p, force


def create_problem(divisions, dt, nonlinear=False):
    model = KM.Model()
    root = model.CreateModelPart('Fluid')
    for variable in (KM.VELOCITY, KM.MESH_VELOCITY, KM.ACCELERATION, KM.PRESSURE,
                     KM.BODY_FORCE, KM.ADVPROJ, KM.DIVPROJ, KM.NODAL_AREA):
        root.AddNodalSolutionStepVariable(variable)
    root.SetBufferSize(3)
    root.ProcessInfo[KM.DOMAIN_SIZE] = 2
    root.ProcessInfo[KM.DYNAMIC_TAU] = 1.0
    root.ProcessInfo[KM.OSS_SWITCH] = 0
    root.ProcessInfo[KM.TIME] = -2*dt
    root.ProcessInfo[KM.DELTA_TIME] = dt
    root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([1.5/dt, -2.0/dt, 0.5/dt])
    properties = root.CreateNewProperties(1)
    properties[KM.DENSITY] = 1.0
    properties[KM.DYNAMIC_VISCOSITY] = MU
    properties[KM.CONSTITUTIVE_LAW] = CFD.Newtonian2DLaw()
    for j in range(divisions+1):
        for i in range(divisions+1):
            node = root.CreateNewNode(j*(divisions+1)+i+1, i/divisions, j/divisions, 0.0)
            for variable in (KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE):
                node.AddDof(variable)
            if i in (0, divisions) or j in (0, divisions):
                node.Fix(KM.VELOCITY_X)
                node.Fix(KM.VELOCITY_Y)
            if i == 0 and j == 0:
                node.Fix(KM.PRESSURE)
    for j in range(divisions):
        for i in range(divisions):
            a = j*(divisions+1)+i+1
            root.CreateNewElement('TimeIntegratedQSVMS2D4N', j*divisions+i+1,
                                  [a, a+1, a+divisions+2, a+divisions+1], properties)
    root.CloneTimeStep(-dt)
    root.CloneTimeStep(0.0)
    for node in root.Nodes:
        for k in range(3):
            u, p, _ = fields(node.X, node.Y, -k*dt, nonlinear)
            node.SetSolutionStepValue(KM.VELOCITY, k, u)
            node.SetSolutionStepValue(KM.PRESSURE, k, p)
    return model, root


def solve(divisions=8, dt=0.02, steps=5, nonlinear=False, bdf1_startup=False):
    model, root = create_problem(divisions, dt, nonlinear)
    # QSVMS does not implement RHS-only assembly: rebuild the combined system
    # each iteration and use an increment criterion, then check the residual.
    strategy = create_qsvms_picard_strategy(root, KM.SkylineLUFactorizationSolver(),
        KM.Parameters('{"relative_tolerance": 1e-9, "absolute_tolerance": 1e-11, '
                      '"maximum_iterations": 50, "echo_level": 0}'))
    bdf = KM.TimeDiscretization.BDF(2)
    strategy.Initialize()
    strategy.Check()
    try:
        for step in range(1, steps+1):
            time = step*dt
            root.CloneTimeStep(time)
            bdf.ComputeAndSaveBDFCoefficients(root.ProcessInfo)
            if bdf1_startup and step == 1:
                root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([1.0/dt, -1.0/dt, 0.0])
            for node in root.Nodes:
                u, p, force = fields(node.X, node.Y, time, nonlinear)
                node.SetSolutionStepValue(KM.BODY_FORCE, force)
                if node.IsFixed(KM.VELOCITY_X):
                    node.SetSolutionStepValue(KM.VELOCITY, u)
                if node.IsFixed(KM.PRESSURE):
                    node.SetSolutionStepValue(KM.PRESSURE, p)
            strategy.InitializeSolutionStep()
            strategy.Predict()
            converged = strategy.SolveSolutionStep()
            strategy.FinalizeSolutionStep()
            if not converged:
                raise RuntimeError('QSVMS nonlinear iteration did not converge')
            residual = {node.Id: [0.0]*3 for node in root.Nodes}
            for element in root.Elements:
                lhs, rhs = KM.Matrix(), KM.Vector()
                element.CalculateLocalSystem(lhs, rhs, root.ProcessInfo)
                for i, node in enumerate(element.GetGeometry()):
                    for d in range(3):
                        residual[node.Id][d] += rhs[3*i+d]
            free_residual = max(abs(residual[node.Id][d]) for node in root.Nodes
                                for d, variable in enumerate((KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE))
                                if not node.IsFixed(variable))
            if free_residual > 1e-8:
                raise RuntimeError(f'QSVMS free-DOF residual is {free_residual}')
        # Integrate the error with the Quad4 two-point tensor Gauss rule.
        u2 = p2 = 0.0
        for element in root.Elements:
            nodes = list(element.GetGeometry())
            hx, hy = nodes[1].X-nodes[0].X, nodes[3].Y-nodes[0].Y
            for xi in (-1/math.sqrt(3), 1/math.sqrt(3)):
                for eta in (-1/math.sqrt(3), 1/math.sqrt(3)):
                    s, t = (1+xi)/2, (1+eta)/2
                    x, y = nodes[0].X+s*hx, nodes[0].Y+t*hy
                    shape = [(1-s)*(1-t), s*(1-t), s*t, (1-s)*t]
                    u = [sum(shape[i]*nodes[i].GetSolutionStepValue(KM.VELOCITY)[d]
                             for i in range(4)) for d in range(2)]
                    p = sum(shape[i]*nodes[i].GetSolutionStepValue(KM.PRESSURE) for i in range(4))
                    ue, pe, _ = fields(x, y, steps*dt, nonlinear)
                    weight = hx*hy/4
                    u2 += weight*sum((u[d]-ue[d])**2 for d in range(2))
                    p2 += weight*(p-pe)**2
        return math.sqrt(u2), math.sqrt(p2)
    finally:
        strategy.Clear()


class TimeIntegratedQSVMSQ1Test(KratosUnittest.TestCase):
    def test_ale_relative_convection_and_pressure_tangent(self):
        _, root = create_problem(2, 0.02)
        root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([75.0, -100.0, 25.0])
        element = next(iter(root.Elements))
        for node in element.GetGeometry():
            for k in range(3):
                node.SetSolutionStepValue(KM.VELOCITY, k, [1+node.X, -node.Y, 0.0])
        element.Initialize(root.ProcessInfo)
        lhs, rhs = KM.Matrix(), KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, root.ProcessInfo)
        # Pressure perturbations leave the frozen convective velocity and tau
        # unchanged, so these columns are true residual derivatives.
        for j, node in enumerate(element.GetGeometry()):
            old = node.GetSolutionStepValue(KM.PRESSURE)
            node.SetSolutionStepValue(KM.PRESSURE, old+1e-6)
            trial_lhs, trial_rhs = KM.Matrix(), KM.Vector()
            element.CalculateLocalSystem(trial_lhs, trial_rhs, root.ProcessInfo)
            node.SetSolutionStepValue(KM.PRESSURE, old)
            for i in range(12):
                self.assertAlmostEqual((rhs[i]-trial_rhs[i])/1e-6, lhs[i, 3*j+2], delta=1e-7)
        # Perturb u and mesh velocity together to keep u-w and tau fixed.
        # This tests the Picard velocity columns, including BDF mass terms.
        for j, node in enumerate(element.GetGeometry()):
            for d in range(2):
                velocity = list(node.GetSolutionStepValue(KM.VELOCITY))
                mesh = list(node.GetSolutionStepValue(KM.MESH_VELOCITY))
                velocity[d] += 1e-6
                mesh[d] += 1e-6
                node.SetSolutionStepValue(KM.VELOCITY, velocity)
                node.SetSolutionStepValue(KM.MESH_VELOCITY, mesh)
                trial_lhs, trial_rhs = KM.Matrix(), KM.Vector()
                element.CalculateLocalSystem(trial_lhs, trial_rhs, root.ProcessInfo)
                velocity[d] -= 1e-6
                mesh[d] -= 1e-6
                node.SetSolutionStepValue(KM.VELOCITY, velocity)
                node.SetSolutionStepValue(KM.MESH_VELOCITY, mesh)
                for i in range(12):
                    self.assertAlmostEqual((rhs[i]-trial_rhs[i])/1e-6, lhs[i, 3*j+d], delta=1e-7)
        # Translating u and w by the same constant preserves relative
        # convection, velocity gradients and BDF acceleration.
        for node in element.GetGeometry():
            for k in range(3):
                u = list(node.GetSolutionStepValue(KM.VELOCITY, k))
                u[0] += 0.2
                node.SetSolutionStepValue(KM.VELOCITY, k, u)
            node.SetSolutionStepValue(KM.MESH_VELOCITY, [0.2, 0.0, 0.0])
        shifted_lhs, shifted_rhs = KM.Matrix(), KM.Vector()
        element.CalculateLocalSystem(shifted_lhs, shifted_rhs, root.ProcessInfo)
        self.assertMatrixAlmostEqual(lhs, shifted_lhs, places=10)
        self.assertVectorAlmostEqual(rhs, shifted_rhs, places=10)

    def test_temporal_refinement(self):
        previous = None
        for dt, steps in ((0.05, 2), (0.025, 4), (0.0125, 8)):
            errors = solve(divisions=8, dt=dt, steps=steps)
            KM.Logger.PrintInfo('QSVMS temporal refinement', f'dt={dt}, errors={errors}')
            if previous:
                self.assertLess(errors[0], previous[0])
                self.assertLess(errors[1], previous[1])
            previous = errors

    def test_spatial_refinement(self):
        previous = None
        for n in (4, 8, 16):
            errors = solve(divisions=n, dt=0.0025, steps=4, nonlinear=True)
            KM.Logger.PrintInfo('QSVMS spatial refinement', f'n={n}, errors={errors}')
            if previous:
                self.assertLess(errors[0], previous[0])
                self.assertLess(errors[1], previous[1])
            previous = errors


if __name__ == '__main__':
    KratosUnittest.main()
