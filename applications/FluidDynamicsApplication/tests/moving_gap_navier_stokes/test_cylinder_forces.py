#    |  /           |
#    ' /   __| _` | __|  _ \   __|
#    . \  |   (   | |   (   |\__ `
#   _|\_\_|  \__,_|\__|\___/ ____/
#                   Multi-Physics
#
#  Fluid Dynamics Application
#  License: BSD License (kratos/license.txt)
#  Main author: Nicolò Antonelli
#
"""Gap-SBM cylinder force tests using independently prescribed fields."""
import math
from types import SimpleNamespace

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.FluidDynamicsApplication as CFD
from KratosMultiphysics.FluidDynamicsApplication.cylinder_sbm_geometry import CylinderSbmGeometry
from KratosMultiphysics.FluidDynamicsApplication.cylinder_traction_utilities import CalculateCylinderPhysicalTraction
from KratosMultiphysics.FluidDynamicsApplication.moving_gap_navier_stokes_solver import MovingGapNavierStokesSolver


RADIUS, GRADIENT, MU, A = .2, 1.0, .3, .4
CENTER = (.5, .5)


def decompose_reaction(boundary, physical):
    """Decompose negative condition RHS; use physical traction independently."""
    penalty = [0.0, 0.0]
    adjoint = [0.0, 0.0]
    for condition in boundary.Conditions:
        nodes = list(condition.GetGeometry())
        a, b, _, d = nodes
        ex, ey, fx, fy = b.X-a.X, b.Y-a.Y, d.X-a.X, d.Y-a.Y
        det = ex*fy-ey*fx
        rule = condition.GetValue(KM.SURROGATE_BOUNDARY_PROJECTION)
        face = condition.GetValue(KM.SURROGATE_BOUNDARY_FACE_COORDINATES)
        prescribed = condition.GetValue(CFD.SBM_BOUNDARY_VELOCITIES)
        alpha = condition.Properties[KM.PENALTY_COEFFICIENT]/math.hypot(
            face[1, 0]-face[0, 0], face[1, 1]-face[0, 1])
        mu = condition.Properties[KM.DYNAMIC_VISCOSITY]
        for q in range(2):
            x, y = rule[q, 0]-a.X, rule[q, 1]-a.Y
            s, t = (x*fy-y*fx)/det, (ex*y-ey*x)/det
            shape = ((1-s)*(1-t), s*(1-t), s*t, (1-s)*t)
            ds, dt = (-(1-t), 1-t, t, -t), (-(1-s), -s, s, 1-s)
            gx = [(ds[i]*fy-dt[i]*ey)/det for i in range(4)]
            gy = [(-ds[i]*fx+dt[i]*ex)/det for i in range(4)]
            velocity = [node.GetSolutionStepValue(KM.VELOCITY) for node in nodes]
            jump = [sum(shape[i]*velocity[i][d] for i in range(4))-prescribed[q, d] for d in range(2)]
            nx, ny, weight = rule[q, 3], rule[q, 4], rule[q, 6]
            for component in range(2):
                penalty[component] += alpha*weight*jump[component]
                for i in range(4):
                    if component == 0:
                        tx = (4*mu*gx[i]/3)*nx + mu*gy[i]*ny
                        ty = mu*gy[i]*nx - (2*mu*gx[i]/3)*ny
                    else:
                        tx = -(2*mu*gy[i]/3)*nx + mu*gx[i]*ny
                        ty = mu*gx[i]*nx + (4*mu*gy[i]/3)*ny
                    adjoint[component] += weight*(tx*jump[0]+ty*jump[1])
    consistency = physical['total'][:2]
    return {'pressure': physical['pressure'][:2], 'viscous': physical['viscous'][:2],
            'consistency': consistency, 'adjoint': adjoint, 'penalty': penalty,
            'total': [consistency[i]+adjoint[i]+penalty[i] for i in range(2)]}


def evaluate(divisions, curvature, pressure_offset=0.0):
    model = KM.Model()
    root = model.CreateModelPart('Fluid')
    for variable in (KM.VELOCITY, KM.PRESSURE, KM.ACCELERATION, KM.BODY_FORCE):
        root.AddNodalSolutionStepVariable(variable)
    root.SetBufferSize(3)
    root.ProcessInfo[KM.DOMAIN_SIZE] = 2
    root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([1., -1., 0.])
    properties = root.CreateNewProperties(1)
    properties[KM.DYNAMIC_VISCOSITY] = MU
    properties[KM.DENSITY] = 1.0
    properties[KM.CONSTITUTIVE_LAW] = CFD.Newtonian2DLaw()
    patch = root.CreateSubModelPart('Patch')
    for j in range(divisions+1):
        for i in range(divisions+1):
            node = patch.CreateNewNode(1+j*(divisions+1)+i, i/divisions, j/divisions, 0.)
            for variable in (KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE):
                node.AddDof(variable)
            node.SetSolutionStepValue(KM.VELOCITY, [curvature*node.Y**2, 0., 0.])
            node.SetSolutionStepValue(KM.PRESSURE, GRADIENT*node.X+pressure_offset)
            # rho=1 and (u.grad)u=0, so f = grad(p)-mu*Laplacian(u).
            node.SetSolutionStepValue(KM.BODY_FORCE, [GRADIENT-2*MU*curvature, 0., 0.])
    for j in range(divisions):
        for i in range(divisions):
            n = 1+j*(divisions+1)+i
            element = patch.CreateNewElement('TimeIntegratedQSVMS2D4N', 1+j*divisions+i,
                                            [n, n+1, n+divisions+2, n+divisions+1], properties)
            element.Set(KM.ACTIVE, True)
    cylinder = CylinderSbmGeometry(patch, CENTER, RADIUS, .05, 10., 'gap_sbm', 256)
    boundary = patch.GetSubModelPart('CylinderBoundary')
    for condition in boundary.Conditions:
        quadrature = condition.GetValue(KM.SURROGATE_BOUNDARY_PROJECTION)
        prescribed = KM.Matrix(2, 3)
        for q in range(2):
            prescribed[q, 0] = curvature*quadrature[q, 1]**2
            prescribed[q, 1] = 0.0
            prescribed[q, 2] = 0.0
        condition.SetValue(CFD.SBM_BOUNDARY_VELOCITIES, prescribed)
        condition.Initialize(root.ProcessInfo)
    # Exercise the solver's public API using these prescribed nodal fields.
    solver = MovingGapNavierStokesSolver.__new__(MovingGapNavierStokesSolver)
    solver._lifecycle = SimpleNamespace(cylinder=cylinder, root=root)
    reaction = solver.CalculateCylinderForce()
    traction = CalculateCylinderPhysicalTraction(boundary)
    decomposition = decompose_reaction(boundary, traction)
    print(f"gap elements={root.GetSubModelPart('CylinderGapVolumes').NumberOfElements()}, "
          f"connectors={root.GetSubModelPart('CylinderGapInterfaces').NumberOfConditions()}, "
          f"dirichlet={boundary.NumberOfConditions()}, decomposition={decomposition}")
    for i in range(2):
        if not math.isclose(reaction[i], decomposition['total'][i], rel_tol=1e-9, abs_tol=1e-10):
            raise AssertionError(f"Condition reaction decomposition mismatch: {reaction}, {decomposition}")
    return reaction, traction


class CylinderForceValidation(KratosUnittest.TestCase):
    def test_pressure_force_and_gauge_invariance(self):
        exact = -GRADIENT*math.pi*RADIUS**2
        errors = []
        for n in (16, 24, 32):
            reaction, traction = evaluate(n, 0.0)
            shifted, shifted_traction = evaluate(n, 0.0, 7.0)
            error = abs(reaction[0]-exact)
            errors.append(error)
            print(f'pressure n={n}: reaction={reaction[:2]}, traction={traction}, '
                  f'abs_error={error:.9g}, relative_error={error/abs(exact):.9g}, '
                  f'reaction_minus_traction={[reaction[i]-traction["total"][i] for i in range(2)]}')
            for i in range(2):
                self.assertAlmostEqual(reaction[i], shifted[i], places=9)
                self.assertAlmostEqual(traction['total'][i], shifted_traction['total'][i], places=9)
            self.assertLess(abs(reaction[1]), .02)
        self.assertLess(errors[-1], errors[0])
        self.assertLess(errors[-1]/abs(exact), .15)

    def test_pressure_and_viscous_traction(self):
        exact = (2*MU*A-GRADIENT)*math.pi*RADIUS**2
        errors = []
        reaction_errors = []
        for n in (16, 24, 32):
            reaction, traction = evaluate(n, A)
            error = abs(traction['total'][0]-exact)
            errors.append(error)
            reaction_errors.append(abs(reaction[0]-exact))
            print(f'viscous n={n}: reaction={reaction[:2]}, pressure={traction["pressure"][:2]}, '
                  f'viscous={traction["viscous"][:2]}, abs_error={error:.9g}, '
                  f'relative_error={error/abs(exact):.9g}, '
                  f'reaction_minus_traction={[reaction[i]-traction["total"][i] for i in range(2)]}')
            self.assertLess(abs(traction['total'][1]), .02)
        self.assertLess(errors[-1], errors[0])
        self.assertLess(errors[-1]/abs(exact), .25)
        self.assertLess(reaction_errors[-1], reaction_errors[0])
        self.assertLess(reaction_errors[-1]/abs(exact), .25)


if __name__ == '__main__':
    KratosUnittest.main()
