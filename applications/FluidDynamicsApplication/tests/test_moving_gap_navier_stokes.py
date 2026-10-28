"""Experimental moving Q1 QSVMS Gap-SBM manufactured solution."""
import json
import math
import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as CFD
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.modelers.multipatch_modeler import FemMultipatchModeler
from KratosMultiphysics.FluidDynamicsApplication.moving_gap_navier_stokes_solver import MovingGapNavierStokesSolver
from test_time_integrated_qsvms_q1 import fields, MU
from test_fluid_patch_gap_coupling import interpolate_nodes
from test_moving_fluid_patch_gap import gradient


def solve_moving(divisions=8, dt=0.01, steps=3, speed=(0.01, 0.0), angle=23.0, penalty=40.0):
    config = {'model_part_name': 'Fluid', 'interface_mode': 'gap_sbm',
              'background': {'origin': [0, 0], 'dimensions': [1, 1],
                             'divisions': [divisions, divisions], 'element_type': 'quadrilateral',
                             'element_name': 'TimeIntegratedQSVMS2D4N'},
              'internal_patches': [{'origin': [0.33, 0.37], 'dimensions': [0.32, 0.32],
                                    'divisions': [max(4, divisions//2)]*2, 'rotation_angle': angle,
                                    'element_type': 'quadrilateral', 'element_name': 'TimeIntegratedQSVMS2D4N'}]}
    model = KM.Model()
    modeler = FemMultipatchModeler(model, KM.Parameters(json.dumps(config)))
    root = modeler.model_part
    for variable in (KM.VELOCITY, KM.MESH_VELOCITY, KM.DISPLACEMENT, KM.ACCELERATION,
                     KM.PRESSURE, KM.BODY_FORCE, KM.ADVPROJ, KM.DIVPROJ, KM.NODAL_AREA):
        root.AddNodalSolutionStepVariable(variable)
    root.SetBufferSize(3)
    root.ProcessInfo[KM.DOMAIN_SIZE] = 2
    root.ProcessInfo[KM.DYNAMIC_TAU] = 1.0
    root.ProcessInfo[KM.OSS_SWITCH] = 0
    root.ProcessInfo[KM.TIME] = 0.0
    root.ProcessInfo[KM.DELTA_TIME] = dt
    root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([1.0/dt, -1.0/dt, 0.0])
    modeler.SetupModelPart()
    for patch in modeler.patches:
        material = next(iter(patch.Elements)).Properties
        material[KM.DENSITY] = 1.0
        material[KM.DYNAMIC_VISCOSITY] = MU
        material[KM.CONSTITUTIVE_LAW] = CFD.Newtonian2DLaw()
    for node in root.Nodes:
        for variable in (KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE):
            node.AddDof(variable)
        velocity, pressure, force = fields(node.X, node.Y, 0.0, True)
        node.SetSolutionStepValue(KM.VELOCITY, velocity)
        node.SetSolutionStepValue(KM.PRESSURE, pressure)
        node.SetSolutionStepValue(KM.BODY_FORCE, force)
        if min(abs(node.X), abs(node.Y), abs(node.X-1), abs(node.Y-1)) < 1e-12:
            node.Fix(KM.VELOCITY_X)
            node.Fix(KM.VELOCITY_Y)
        if abs(node.X) < 1e-12 and abs(node.Y) < 1e-12:
            node.Fix(KM.PRESSURE)
    # Establish previous ProcessInfo for the BDF utility. Startup uses BDF1,
    # so no negative-time numerical velocity is consumed.
    root.ProcessInfo[KM.TIME] = -dt
    root.CloneTimeStep(0.0)
    props = root.CreateNewProperties(100)
    props[KM.PENALTY_COEFFICIENT] = penalty
    coupling = {'background_model_part_name': 'Fluid.Patch_0', 'inner_model_part_name': 'Fluid.Patch_1',
                'interface_model_part_name': 'Fluid.Coupling', 'properties_id': 100,
                'allow_transient': True, 'formulation': 'qsvms'}
    CFD.CreateFluidPatchGapCouplingConditionsProcess(model, KM.Parameters(json.dumps(coupling))).Execute()
    computing = root.CreateSubModelPart('Computing')
    computing.AddElements([e.Id for e in root.Elements if e.Is(KM.ACTIVE)])
    computing.AddConditions([c.Id for c in root.GetSubModelPart('Coupling').Conditions])
    ids = {n.Id for e in computing.Elements for n in e.GetGeometry()}
    ids.update(n.Id for c in computing.Conditions for n in c.GetGeometry())
    computing.AddNodes(sorted(ids))
    settings = {'time_step': dt, 'translation_velocity': list(speed), 'bdf1_startup': True,
                'coupling_settings': coupling, 'relative_tolerance': 1e-9,
                'absolute_tolerance': 1e-11, 'maximum_iterations': 50}
    solver = MovingGapNavierStokesSolver(model, KM.Parameters(json.dumps(settings)))
    solver.Initialize()
    rows = []
    for step in range(steps):
        time = solver.AdvanceInTime(step*dt)
        computing = solver.GetComputingModelPart()
        for node in computing.Nodes:
            u, p, force = fields(node.X, node.Y, time, True)
            node.SetSolutionStepValue(KM.BODY_FORCE, force)
            if node.IsFixed(KM.VELOCITY_X):
                node.SetSolutionStepValue(KM.VELOCITY, u)
            if node.IsFixed(KM.PRESSURE):
                node.SetSolutionStepValue(KM.PRESSURE, p)
        solver.InitializeSolutionStep()
        if not solver.SolveSolutionStep():
            raise RuntimeError('Moving QSVMS Picard iteration failed to converge')
        l2u = l2p = area = div2 = signed_divergence = 0.0
        for element_id, point, weight in KM.FemMultipatchUtilities.GetGaussPointsWithWeights(computing):
            u, p = interpolate_nodes(list(computing.GetElement(element_id).GetGeometry()), point)
            ue, pe, _ = fields(point[0], point[1], time, True)
            l2u += weight*sum((u[d]-ue[d])**2 for d in range(2))
            l2p += weight*(p-pe)**2
            area += weight
            grad = gradient(list(computing.GetElement(element_id).GetGeometry()), point)
            divergence = grad[0][0]+grad[1][1]
            div2 += weight*divergence*divergence
            signed_divergence += weight*divergence
        jump2 = 0.0
        for condition in root.GetSubModelPart('Coupling.PhysicalInterfaces').Conditions:
            nodes = list(condition.GetGeometry())
            points = condition.CalculateOnIntegrationPoints(KM.INTEGRATION_COORDINATES, root.ProcessInfo)
            length = math.sqrt(3)*math.dist(points[0], points[1])
            for point in points:
                ua, _ = interpolate_nodes(nodes[:4], point)
                ub, _ = interpolate_nodes(nodes[4:], point)
                jump2 += 0.5*length*sum((ua[d]-ub[d])**2 for d in range(2))
        pressure_max_error = max(abs(node.GetSolutionStepValue(KM.PRESSURE)-
                                     fields(node.X, node.Y, time, True)[1]) for node in computing.Nodes)
        row = {'time': time, 'velocity_l2': math.sqrt(l2u), 'pressure_l2': math.sqrt(l2p),
               'area': area, 'divergence_l2': math.sqrt(div2),
               'integrated_divergence': signed_divergence, 'interface_jump_l2': math.sqrt(jump2),
               'pressure_max_error': pressure_max_error,
               'free_dof_residual': solver.free_dof_residual, **solver.motion_diagnostics}
        KM.Logger.PrintInfo('MovingGapQSVMS', str(row))
        rows.append(row)
        solver.FinalizeSolutionStep()
    return rows


class MovingGapNavierStokesTest(KratosUnittest.TestCase):
    def test_moving_manufactured_and_refinement(self):
        previous = None
        for n in (8, 12, 16):
            result = solve_moving(divisions=n)[-1]
            self.assertAlmostEqual(result['area'], 1.0, delta=1e-9)
            self.assertTrue(math.isfinite(result['velocity_l2']) and math.isfinite(result['pressure_l2']))
            self.assertLess(result['free_dof_residual'], 1e-8)
            if previous is not None:
                self.assertLess(result['velocity_l2'], previous['velocity_l2'])
                self.assertLess(result['pressure_l2'], previous['pressure_l2'])
            previous = result


if __name__ == '__main__':
    KratosUnittest.main()
