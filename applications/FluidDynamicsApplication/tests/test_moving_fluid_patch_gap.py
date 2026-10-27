"""Transient Stokes MMS: translating ALE square and fixed Eulerian background.

Run this file with KratosUnittest, or import solve_moving_gap for diagnostics.
Reactivation histories can be transferred from the inner patch or seeded from
the exact field for comparison. No interface solution values are prescribed.
"""
import math
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication.rebuild_fluid_patch_gap_process import RebuildFluidPatchGapProcess
from test_fluid_patch_gap_coupling import create_gap_problem, exact_fields, interpolate_nodes


def transient_fields(x, y, time, nonlinear=False):
    velocity, pressure, steady_force = exact_fields(x, y, nonlinear)
    amplitude = 1.0 + time
    # Eulerian forcing: du/dt - mu*Laplacian(u) + grad(p), rho=1.
    # No mesh transport is put in f: the ALE derivative supplies +w.grad(u).
    return ([amplitude*v for v in velocity], amplitude*pressure,
            [velocity[d] + amplitude*steady_force[d] for d in range(3)])


def seed_history(node, time, dt, speed, nonlinear):
    for step in range(3):
        u, p, _ = transient_fields(node.X-speed*step*dt, node.Y, time-step*dt, nonlinear)
        node.SetSolutionStepValue(KM.VELOCITY, step, u)
        node.SetSolutionStepValue(KM.PRESSURE, step, p)


def transfer_background_histories(inner, nodes, time, dt, speed, nonlinear, fallback=False):
    """Transfer steps 1 and 2 from historical ALE positions, without extrapolation.

    The locator sees the current mesh. For rigid translation, translating the
    fixed query forward by speed*k*dt gives the same Q1 local coordinates as
    searching the actual mesh at t-k*dt. Collect all results before writing.
    """
    locator = KM.BinBasedFastPointLocator2D(inner)
    locator.UpdateSearchDatabase()
    values, missing = {}, []
    for node in nodes:
        for k in (1, 2):
            found, shape, element = locator.FindPointOnMesh([node.X + speed*k*dt, node.Y, 0.0])
            if found:
                source = list(element.GetGeometry())
                values[node.Id, k] = [sum(shape[i]*source[i].GetSolutionStepValue(KM.VELOCITY, k)[d]
                                           for i in range(4)) for d in range(3)]
            else:
                missing.append((node.Id, k, time-k*dt, node.X, node.Y))
                if fallback:
                    values[node.Id, k] = transient_fields(node.X, node.Y, time-k*dt, nonlinear)[0]
    if missing:
        message = f"Historical inner patch does not contain (node, step, time, x, y): {missing}"
        if not fallback:
            raise RuntimeError(message)
        KM.Logger.PrintWarning("MovingGapTransfer", message + "; using explicitly requested exact fallback.")
    for node in nodes:
        for k in (1, 2):
            node.SetSolutionStepValue(KM.VELOCITY, k, values[node.Id, k])
        # A numerical starting guess; BDF2 uses slots 1 and 2 as its histories.
        node.SetSolutionStepValue(KM.VELOCITY, values[node.Id, 1])
    return len(missing)


def gradient(nodes, point):
    """Physical Q1 gradient on the affine quad, including gap extrapolation."""
    a, b, d = nodes[0], nodes[1], nodes[3]
    ex, ey, fx, fy = b.X-a.X, b.Y-a.Y, d.X-a.X, d.Y-a.Y
    determinant = ex*fy-ey*fx
    s = ((point[0]-a.X)*fy-(point[1]-a.Y)*fx)/determinant
    t = (ex*(point[1]-a.Y)-ey*(point[0]-a.X))/determinant
    ds, dt = [t-1, 1-t, t, -t], [s-1, -s, s, 1-s]
    result = [[0.0, 0.0], [0.0, 0.0]]
    for i, node in enumerate(nodes):
        gx, gy = (fy*ds[i]-ey*dt[i])/determinant, (-fx*ds[i]+ex*dt[i])/determinant
        for component in range(2):
            value = node.GetSolutionStepValue(KM.VELOCITY)[component]
            result[component][0] += gx*value
            result[component][1] += gy*value
    return result


def diagnostics(root, computing, time, nonlinear):
    l2u = l2p = div2 = volume = signed_div = 0.0
    max_pressure_error = 0.0
    for node in computing.Nodes:
        _, p, _ = transient_fields(node.X, node.Y, time, nonlinear)
        max_pressure_error = max(max_pressure_error, abs(node.GetSolutionStepValue(KM.PRESSURE)-p))
    for element_id, point, weight in KM.FemMultipatchUtilities.GetGaussPointsWithWeights(computing):
        nodes = list(computing.GetElement(element_id).GetGeometry())
        u, p = interpolate_nodes(nodes, point)
        ue, pe, _ = transient_fields(point[0], point[1], time, nonlinear)
        grad = gradient(nodes, point)
        div = grad[0][0]+grad[1][1]
        l2u += weight*sum((u[d]-ue[d])**2 for d in range(2))
        l2p += weight*(p-pe)**2
        div2 += weight*div*div
        signed_div += weight*div
        volume += weight
    # Broken divergence alone misses flux jumps at interfaces.
    jump2 = 0.0
    physical = root.GetSubModelPart("Coupling.PhysicalInterfaces")
    for condition in physical.Conditions:
        nodes = list(condition.GetGeometry())
        points = condition.CalculateOnIntegrationPoints(KM.INTEGRATION_COORDINATES, root.ProcessInfo)
        # Two equally weighted Gauss points on this physical segment.
        length = math.sqrt(3)*math.dist(points[0], points[1])
        for point in points:
            ua, _ = interpolate_nodes(nodes[:4], point)
            ub, _ = interpolate_nodes(nodes[4:], point)
            jump2 += 0.5*length*sum((ua[d]-ub[d])**2 for d in range(2))
    return dict(velocity_l2=math.sqrt(l2u), pressure_l2=math.sqrt(l2p),
                divergence_l2=math.sqrt(div2), integrated_divergence=signed_div,
                interface_jump_l2=math.sqrt(jump2), pressure_spike=max_pressure_error,
                integrated_area=volume)


def solve_moving_gap(divisions=12, dt=0.02, steps=6, speed=0.01,
                     coefficient=0.0, nonlinear=False, angle=0.0,
                     history_mode="exact"):
    if history_mode not in ("exact", "transfer", "transfer_with_exact_fallback"):
        raise ValueError("Unknown history mode: " + history_mode)
    model, modeler, computing, settings = create_gap_problem(
        divisions=divisions, inner_divisions=max(4, divisions//2),
        coefficient=coefficient, nonlinear=nonlinear, angle=angle, ale=True,
        patch_origin=[0.249, 0.37])
    root = modeler.model_part
    background, inner = modeler.patches
    inner_ids = {n.Id for n in inner.Nodes}
    connectivity = {e.Id: tuple(n.Id for n in e.GetGeometry()) for e in inner.Elements}
    node_count = root.NumberOfNodes()
    root.ProcessInfo[KM.TIME] = -2*dt
    root.CloneTimeStep(-dt)
    root.CloneTimeStep(0.0)
    root.ProcessInfo[KM.DYNAMIC_TAU] = 0.0  # Same quasi-static ASGS as the steady MMS.
    for node in root.Nodes:
        seed_history(node, 0.0, dt, speed if node.Id in inner_ids else 0.0, nonlinear)
        local_speed = speed if node.Id in inner_ids else 0.0
        for k in range(3):
            node.SetSolutionStepValue(KM.DISPLACEMENT, k, [-local_speed*k*dt, 0.0, 0.0])
            node.SetSolutionStepValue(KM.MESH_VELOCITY, k, [local_speed, 0.0, 0.0])
    bdf = KM.TimeDiscretization.BDF(2)
    rebuild = RebuildFluidPatchGapProcess(model, [settings])
    previous_nodes = {n.Id for n in computing.Nodes}
    previous_elements = {e.Id for e in background.Elements if e.Is(KM.ACTIVE)}
    previous_triangles = root.GetSubModelPart("Interfaces.Patch_1.GapTriangles").NumberOfGeometries()
    rows = []
    for step in range(1, steps+1):
        time = step*dt
        root.CloneTimeStep(time)
        bdf.ComputeAndSaveBDFCoefficients(root.ProcessInfo)
        for node in inner.Nodes:
            node.X = node.X0 + speed*time
            node.Y = node.Y0
            node.SetSolutionStepValue(KM.DISPLACEMENT, [speed*time, 0.0, 0.0])
            # Same BDF operator for mesh displacement and fluid velocity.
            mesh_speed = sum(root.ProcessInfo[KM.BDF_COEFFICIENTS][k]*
                             node.GetSolutionStepValue(KM.DISPLACEMENT, k)[0] for k in range(3))
            node.SetSolutionStepValue(KM.MESH_VELOCITY, [mesh_speed, 0.0, 0.0])
        rebuild.Execute()
        computing = rebuild.computing
        current_nodes = {n.Id for n in computing.Nodes}
        current_elements = {e.Id for e in background.Elements if e.Is(KM.ACTIVE)}
        new_nodes = current_nodes-previous_nodes-inner_ids
        # Histories of already-active nodes are untouched after CloneTimeStep.
        retained = {node_id: [tuple(root.GetNode(node_id).GetSolutionStepValue(KM.VELOCITY, k))
                              for k in (1, 2)] for node_id in current_nodes & previous_nodes}
        targets = [root.GetNode(node_id) for node_id in sorted(new_nodes)]
        missing = 0
        if history_mode == "exact":
            for node in targets:
                seed_history(node, time, dt, 0.0, nonlinear)
        elif targets:
            missing = transfer_background_histories(inner, targets, time, dt, speed, nonlinear,
                fallback=history_mode == "transfer_with_exact_fallback")
        if any(retained[node_id] != [tuple(root.GetNode(node_id).GetSolutionStepValue(KM.VELOCITY, k))
                                    for k in (1, 2)] for node_id in retained):
            raise RuntimeError("Reactivation overwrote a valid BDF history.")
        for node in computing.Nodes:
            u, p, force = transient_fields(node.X, node.Y, time, nonlinear)
            node.SetSolutionStepValue(KM.BODY_FORCE, force)
            if node.IsFixed(KM.VELOCITY_X):
                node.SetSolutionStepValue(KM.VELOCITY, u)
            if node.IsFixed(KM.PRESSURE):
                node.SetSolutionStepValue(KM.PRESSURE, p)
        strategy = KM.ResidualBasedLinearStrategy(
            computing, KM.ResidualBasedIncrementalUpdateStaticScheme(),
            KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver()),
            False, False, False, False)
        strategy.SetEchoLevel(0)
        try:
            strategy.Initialize()
            strategy.Check()
            strategy.Solve()
            row = diagnostics(root, computing, time, nonlinear)
            # Assemble a free-DOF residual independently of the strategy.
            residual = {n.Id: [0.0]*3 for n in computing.Nodes}
            for entity in [*computing.Elements, *computing.Conditions]:
                lhs, rhs = KM.Matrix(), KM.Vector()
                entity.CalculateLocalSystem(lhs, rhs, root.ProcessInfo)
                for i, node in enumerate(entity.GetGeometry()):
                    for d in range(3):
                        residual[node.Id][d] += rhs[3*i+d]
            row["free_residual"] = max(abs(residual[n.Id][d]) for n in computing.Nodes
                for d, variable in enumerate((KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE)) if not n.IsFixed(variable))
        finally:
            # No strategy may retain references to entities deleted next step.
            strategy.Clear()
        triangles = root.GetSubModelPart("Interfaces.Patch_1.GapTriangles")
        row.update(time=time, newly_active_nodes=len(new_nodes),
                   missing_history_samples=missing,
                   transferred_history_samples=0 if history_mode == "exact" else 2*len(new_nodes)-missing,
                   newly_active_elements=len(current_elements-previous_elements),
                   deactivated_elements=len(previous_elements-current_elements),
                   gap_triangles=triangles.NumberOfGeometries(),
                   gap_triangle_change=triangles.NumberOfGeometries()-previous_triangles,
                   gap_area=sum(g.Area() for g in triangles.Geometries))
        if root.NumberOfNodes() != node_count:
            raise RuntimeError("Motion introduced solution nodes.")
        if connectivity != {e.Id: tuple(n.Id for n in e.GetGeometry()) for e in inner.Elements}:
            raise RuntimeError("Rigid motion changed inner connectivity.")
        if any(n.X != n.X0 or n.Y != n.Y0 for n in background.Nodes):
            raise RuntimeError("Eulerian background moved.")
        KM.Logger.PrintInfo("MovingGapStokes", str(row))
        rows.append(row)
        previous_nodes, previous_elements = current_nodes, current_elements
        previous_triangles = triangles.NumberOfGeometries()
    return rows


class MovingFluidPatchGapTest(KratosUnittest.TestCase):
    def test_historical_q1_transfer_and_outside_rejection(self):
        _, modeler, _, _ = create_gap_problem(ale=True, angle=23.0)
        root = modeler.model_part
        inner = modeler.patches[1]
        speed, dt, time = 0.01, 0.02, 0.12
        for node in inner.Nodes:
            for k in (1, 2):
                x, y = node.X-speed*k*dt, node.Y
                node.SetSolutionStepValue(KM.VELOCITY, k, [2*x+3*y+k, -x+4*y-k, 0.0])
        parent = next(iter(inner.Elements))
        center = parent.GetGeometry().Center()
        new_id = max(n.Id for n in root.Nodes)+1
        target = root.CreateNewNode(new_id, center.X-0.5*speed*dt, center.Y, 0.0)
        self.assertEqual(transfer_background_histories(inner, [target], time, dt, speed, False), 0)
        for k in (1, 2):
            u = target.GetSolutionStepValue(KM.VELOCITY, k)
            self.assertAlmostEqual(u[0], 2*target.X+3*target.Y+k, delta=1e-10)
            self.assertAlmostEqual(u[1], -target.X+4*target.Y-k, delta=1e-10)
        outside = root.CreateNewNode(new_id+1, -1.0, -1.0, 0.0)
        target.SetSolutionStepValue(KM.VELOCITY, 1, [9.0, 8.0, 0.0])
        with self.assertRaisesRegex(RuntimeError, "Historical inner patch does not contain"):
            transfer_background_histories(inner, [target, outside], time, dt, speed, False)
        self.assertAlmostEqual(target.GetSolutionStepValue(KM.VELOCITY, 1)[0], 9.0)

    def test_compare_transferred_and_exact_histories(self):
        for nonlinear in (False, True):
            exact = solve_moving_gap(nonlinear=nonlinear, history_mode="exact")
            numerical = solve_moving_gap(nonlinear=nonlinear, history_mode="transfer_with_exact_fallback")
            self.assertEqual([r["newly_active_nodes"] for r in exact],
                             [r["newly_active_nodes"] for r in numerical])
            self.assertGreater(sum(r["newly_active_nodes"] for r in numerical), 0)
            self.assertGreater(sum(r["transferred_history_samples"] for r in numerical), 0)
            for a, b in zip(exact, numerical):
                self.assertTrue(all(math.isfinite(b[key]) for key in
                                    ("velocity_l2", "pressure_l2", "pressure_spike", "free_residual")))
                self.assertLess(b["free_residual"], 1e-8)
                KM.Logger.PrintInfo("MovingGapHistoryComparison", f"exact={a}; numerical={b}")

    def test_ale_local_linearization_and_zero_mesh_velocity(self):
        _, modeler, _, _ = create_gap_problem(ale=True)
        root = modeler.model_part
        element = next(iter(modeler.patches[1].Elements))
        nodes = list(element.GetGeometry())
        root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([75.0, -100.0, 25.0])
        root.ProcessInfo[KM.DELTA_TIME] = 0.02
        root.ProcessInfo[KM.DYNAMIC_TAU] = 1.0
        for node in nodes:
            seed_history(node, 0.0, 0.02, 0.01, False)
            node.SetSolutionStepValue(KM.MESH_VELOCITY, [0.01, 0.0, 0.0])
        element.Initialize(root.ProcessInfo)
        lhs, rhs = KM.Matrix(), KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, root.ProcessInfo)
        epsilon = 1e-6
        for j, node in enumerate(nodes):
            for d in range(2):
                value = list(node.GetSolutionStepValue(KM.VELOCITY))
                value[d] += epsilon
                node.SetSolutionStepValue(KM.VELOCITY, value)
                perturbed = KM.Vector()
                element.CalculateRightHandSide(perturbed, root.ProcessInfo)
                value[d] -= epsilon
                node.SetSolutionStepValue(KM.VELOCITY, value)
                for i in range(12):
                    self.assertAlmostEqual((rhs[i]-perturbed[i])/epsilon, lhs[i, 3*j+d], delta=1e-7)
        for node in nodes:
            node.SetSolutionStepValue(KM.MESH_VELOCITY, [0.0, 0.0, 0.0])
        twin = root.CreateNewElement("SymbolicStokes2D4N", max(e.Id for e in root.Elements)+1,
                                     [n.Id for n in nodes], element.Properties)
        twin.Initialize(root.ProcessInfo)
        base_lhs, base_rhs = KM.Matrix(), KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, root.ProcessInfo)
        twin.CalculateLocalSystem(base_lhs, base_rhs, root.ProcessInfo)
        self.assertMatrixAlmostEqual(lhs, base_lhs, places=12)
        self.assertVectorAlmostEqual(rhs, base_rhs, places=12)

    def test_linear_moving_and_stationary(self):
        for speed in (0.0, 0.01):
            for coefficient in (0.0, 40.0):
                rows = solve_moving_gap(speed=speed, coefficient=coefficient)
                for row in rows:
                    for name in ("velocity_l2", "pressure_l2", "interface_jump_l2", "pressure_spike",
                                 "divergence_l2", "free_residual"):
                        self.assertLess(row[name], 1e-8, str(row))
                    self.assertAlmostEqual(row["integrated_area"], 1.0, delta=1e-10)
                if speed:
                    self.assertGreater(sum(r["newly_active_elements"] for r in rows), 0)

    def test_nonlinear_spatial_refinement(self):
        for coefficient in (0.0, 40.0):
            previous = None
            for divisions in (8, 16, 32):
                # Refine time too, keeping T=0.12 fixed, to avoid a BDF2 floor.
                steps = 6*(divisions//8)
                row = solve_moving_gap(divisions=divisions, dt=0.12/steps, steps=steps,
                                       nonlinear=True, coefficient=coefficient)[-1]
                self.assertTrue(all(math.isfinite(value) for value in row.values()))
                self.assertLess(row["free_residual"], 1e-8)
                self.assertAlmostEqual(row["integrated_area"], 1.0, delta=1e-10)
                if previous:
                    self.assertLess(row["velocity_l2"], previous["velocity_l2"])
                    self.assertLess(row["pressure_l2"], previous["pressure_l2"])
                previous = row


if __name__ == "__main__":
    KratosUnittest.main()
