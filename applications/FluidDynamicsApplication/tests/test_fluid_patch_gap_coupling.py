"""Assembled steady Stokes tests for positive-gap Q1 multipatch coupling.

No exact values are prescribed on either side of the internal interface.
Only exterior velocities and one global pressure reference are fixed.
"""
import json
import math
import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as CFD
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.modelers.multipatch_modeler import FemMultipatchModeler


def exact_fields(x, y, nonlinear=False):
    if nonlinear:
        # Stream function x^2*y^2 gives a divergence-free cubic velocity.
        return [2*x*x*y, -2*x*y*y, 0.0], x*x+y*y, [2*x-4*1.3*y, 2*y+4*1.3*x, 0.0]
    return [1+2*x+3*y, -2+4*x-2*y, 0.0], 1+2*x-3*y, [2.0, -3.0, 0.0]


def interpolate(element, point):
    """Independent affine Q1 interpolation, used only for test diagnostics."""
    return interpolate_nodes(list(element.GetGeometry()), point)


def interpolate_nodes(nodes, point):
    a, b, d = nodes[0], nodes[1], nodes[3]
    ex, ey = b.X-a.X, b.Y-a.Y
    fx, fy = d.X-a.X, d.Y-a.Y
    dx, dy = point[0]-a.X, point[1]-a.Y
    determinant = ex*fy-ey*fx
    s, t = (dx*fy-dy*fx)/determinant, (ex*dy-ey*dx)/determinant
    shape = [(1-s)*(1-t), s*(1-t), s*t, (1-s)*t]
    velocity = [sum(shape[i]*nodes[i].GetSolutionStepValue(KM.VELOCITY)[k]
                    for i in range(4)) for k in range(2)]
    pressure = sum(shape[i]*nodes[i].GetSolutionStepValue(KM.PRESSURE) for i in range(4))
    return velocity, pressure


def create_gap_problem(angle=23.0, coefficient=0.0, divisions=12,
                       inner_divisions=7, nonlinear=False, fitted=False, ale=False,
                       patch_origin=None):
    origin, size = ([0.25, 0.25], 0.5) if fitted else ([0.33, 0.37], 0.32)
    if patch_origin is not None:
        origin = patch_origin
    settings = KM.Parameters(json.dumps({
        "model_part_name": "Fluid", "interface_mode": "gap_sbm",
        "background": {"origin": [0, 0], "dimensions": [1, 1],
                       "divisions": [divisions, divisions], "element_type": "quadrilateral",
                       "element_name": "SymbolicStokes2D4N"},
        "internal_patches": [{"origin": origin, "dimensions": [size, size],
                              "divisions": [inner_divisions, inner_divisions+1],
                              "element_type": "quadrilateral", "rotation_angle": angle,
                              "element_name": "SymbolicStokesALE2D4N" if ale else "SymbolicStokes2D4N"}]}))
    model = KM.Model()
    modeler = FemMultipatchModeler(model, settings)
    root = modeler.model_part
    for variable in (KM.VELOCITY, KM.ACCELERATION, KM.PRESSURE, KM.BODY_FORCE):
        root.AddNodalSolutionStepVariable(variable)
    if ale:
        root.AddNodalSolutionStepVariable(KM.MESH_VELOCITY)
        root.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    root.SetBufferSize(3)
    root.ProcessInfo[KM.DELTA_TIME] = 1.0
    root.ProcessInfo[KM.DYNAMIC_TAU] = 0.0
    root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([0.0, 0.0, 0.0])
    modeler.SetupModelPart()
    for patch in modeler.patches:
        material = next(iter(patch.Elements)).Properties
        material[KM.DENSITY] = 1.0
        material[KM.DYNAMIC_VISCOSITY] = 1.3
        material[KM.CONSTITUTIVE_LAW] = CFD.Newtonian2DLaw()
    for node in root.Nodes:
        for variable in (KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE):
            node.AddDof(variable)
        velocity, pressure, force = exact_fields(node.X, node.Y, nonlinear)
        node.SetSolutionStepValue(KM.BODY_FORCE, force)
        if min(abs(node.X), abs(node.Y), abs(node.X-1), abs(node.Y-1)) < 1e-12:
            node.SetSolutionStepValue(KM.VELOCITY, velocity)
            node.Fix(KM.VELOCITY_X)
            node.Fix(KM.VELOCITY_Y)
        if abs(node.X) < 1e-12 and abs(node.Y) < 1e-12:
            node.SetSolutionStepValue(KM.PRESSURE, pressure)
            node.Fix(KM.PRESSURE)
    properties = root.CreateNewProperties(100)
    properties[KM.PENALTY_COEFFICIENT] = coefficient
    process_settings = KM.Parameters(json.dumps({
        "background_model_part_name": "Fluid.Patch_0", "inner_model_part_name": "Fluid.Patch_1",
        "interface_model_part_name": "Fluid.Coupling", "properties_id": 100}))
    node_count = root.NumberOfNodes()
    CFD.CreateFluidPatchGapCouplingConditionsProcess(model, process_settings).Execute()
    if root.NumberOfNodes() != node_count:
        raise RuntimeError("Gap coupling introduced solution nodes.")
    computing = root.CreateSubModelPart("Computing")
    computing.AddElements([e.Id for e in root.Elements if e.Is(KM.ACTIVE)])
    computing.AddConditions([c.Id for c in root.GetSubModelPart("Coupling").Conditions])
    ids = {n.Id for e in computing.Elements for n in e.GetGeometry()}
    ids.update(n.Id for c in computing.Conditions for n in c.GetGeometry())
    computing.AddNodes(sorted(ids))
    return model, modeler, computing, process_settings


def solve_gap_problem(**kwargs):
    model, modeler, computing, settings = create_gap_problem(**kwargs)
    strategy = KM.ResidualBasedLinearStrategy(
        computing, KM.ResidualBasedIncrementalUpdateStaticScheme(),
        KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver()),
        False, False, False, False)
    strategy.SetEchoLevel(0)
    strategy.Initialize()
    strategy.Check()
    strategy.Solve()
    nonlinear = kwargs.get("nonlinear", False)
    errors = {"velocity": 0.0, "pressure": 0.0, "velocity_jump": 0.0, "pressure_jump": 0.0}
    for node in computing.Nodes:
        u, p, _ = exact_fields(node.X, node.Y, nonlinear)
        errors["velocity"] = max(errors["velocity"], *(abs(node.GetSolutionStepValue(KM.VELOCITY)[d]-u[d]) for d in range(2)))
        errors["pressure"] = max(errors["pressure"], abs(node.GetSolutionStepValue(KM.PRESSURE)-p))
    for c in modeler.model_part.GetSubModelPart("Coupling.PhysicalInterfaces").Conditions:
        nodes = list(c.GetGeometry())
        for point in c.CalculateOnIntegrationPoints(KM.INTEGRATION_COORDINATES, computing.ProcessInfo):
            a, b = interpolate_nodes(nodes[:4], point), interpolate_nodes(nodes[4:], point)
            errors["velocity_jump"] = max(errors["velocity_jump"], *(abs(a[0][d]-b[0][d]) for d in range(2)))
            errors["pressure_jump"] = max(errors["pressure_jump"], abs(a[1]-b[1]))
    l2u = l2p = 0.0
    for element_id, point, weight in KM.FemMultipatchUtilities.GetGaussPointsWithWeights(computing):
        u, p = interpolate(computing.GetElement(element_id), point)
        ue, pe, _ = exact_fields(point[0], point[1], nonlinear)
        l2u += weight*sum((u[d]-ue[d])**2 for d in range(2))
        l2p += weight*(p-pe)**2
    errors.update(velocity_l2=math.sqrt(l2u), pressure_l2=math.sqrt(l2p))
    strategy.Clear()
    return model, modeler, computing, errors


class FluidPatchGapCouplingTest(KratosUnittest.TestCase):
    def test_linear_translated_and_rotated(self):
        for angle in (0.0, 23.0, 45.0):
            for coefficient in (0.0, 40.0):
                with self.subTest(angle=angle, coefficient=coefficient):
                    _, modeler, computing, errors = solve_gap_problem(angle=angle, coefficient=coefficient)
                    for name, value in errors.items():
                        self.assertTrue(math.isfinite(value), name)
                        self.assertLess(value, 1e-8, f"{name}: {errors}")
                    self.assertEqual(sum(n.IsFixed(KM.PRESSURE) for n in computing.Nodes), 1)
                    for n in modeler.patches[1].Nodes:
                        self.assertFalse(n.IsFixed(KM.VELOCITY_X))
                        self.assertFalse(n.IsFixed(KM.PRESSURE))
                    KM.Logger.PrintInfo("FluidPatchGapLinear", f"angle={angle}, coefficient={coefficient}, {errors}")

    def test_fitted_zero_gap_limit(self):
        _, modeler, _, errors = solve_gap_problem(angle=0.0, fitted=True)
        self.assertEqual(modeler.model_part.GetSubModelPart("Coupling.GapVolumes").NumberOfElements(), 0)
        self.assertLess(max(errors.values()), 1e-8)

    def test_exact_residual_balance(self):
        for coefficient in (0.0, 40.0):
            model, modeler, computing, _ = create_gap_problem(coefficient=coefficient)
            for node in computing.Nodes:
                u, p, _ = exact_fields(node.X, node.Y)
                node.SetSolutionStepValue(KM.VELOCITY, u)
                node.SetSolutionStepValue(KM.PRESSURE, p)
            residual = {n.Id: [0.0]*3 for n in computing.Nodes}
            for entity in list(computing.Elements) + list(computing.Conditions):
                entity.Initialize(computing.ProcessInfo)
                lhs, rhs = KM.Matrix(), KM.Vector()
                entity.CalculateLocalSystem(lhs, rhs, computing.ProcessInfo)
                for i, node in enumerate(entity.GetGeometry()):
                    for d in range(3):
                        residual[node.Id][d] += rhs[3*i+d]
            for node in computing.Nodes:
                for d, variable in enumerate((KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE)):
                    if not node.IsFixed(variable):
                        self.assertAlmostEqual(residual[node.Id][d], 0.0, delta=1e-10)

    def test_creator_rejects_duplicate_execution(self):
        model, _, _, settings = create_gap_problem()
        with self.assertRaisesRegex(RuntimeError, "only once"):
            CFD.CreateFluidPatchGapCouplingConditionsProcess(model, settings).Execute()

    def test_nonlinear_manufactured_refinement(self):
        for coefficient in (0.0, 40.0):
            previous = None
            for divisions, inner in ((8, 4), (16, 8), (32, 16)):
                _, _, _, errors = solve_gap_problem(divisions=divisions, inner_divisions=inner,
                                                    coefficient=coefficient, nonlinear=True)
                if previous:
                    self.assertLess(errors["velocity_l2"], previous["velocity_l2"])
                    self.assertLess(errors["pressure_l2"], previous["pressure_l2"])
                KM.Logger.PrintInfo("FluidPatchGapRefinement", f"n={divisions}, coefficient={coefficient}, {errors}")
                previous = errors


if __name__ == "__main__":
    KratosUnittest.main()
