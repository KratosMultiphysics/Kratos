#   KRATOS  _____________
#          /  _/ ____/   |
#          / // / __/ /| |
#        _/ // /_/ / ___ |
#       /___/\____/_/  |_| Application

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import unittest


class Shell6pElementTests(KratosUnittest.TestCase):

    def solve_cantilever(create_geometry):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

        shell_properties = model_part.GetProperties()[1]
        shell_properties.SetValue(KM.THICKNESS, 0.1)
        shell_properties.SetValue(KM.YOUNG_MODULUS, 200000000)
        shell_properties.SetValue(KM.POISSON_RATIO, 0)
        shell_properties.SetValue(KM.CONSTITUTIVE_LAW, SMA.LinearElastic3DLaw())

        surface = create_geometry(model_part)

        quadrature_point_geometries = KM.GeometriesVector()

        surface.CreateQuadraturePointGeometries(quadrature_point_geometries, 3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('Shell6pElement', element_id, quadrature_point_geometries[i], shell_properties)
            element_id += 1

        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)

        surface[0].Fix(KM.DISPLACEMENT_X)
        surface[0].Fix(KM.DISPLACEMENT_Y)
        surface[0].Fix(KM.DISPLACEMENT_Z)

        surface[3].Fix(KM.DISPLACEMENT_X)
        surface[3].Fix(KM.DISPLACEMENT_Y)
        surface[3].Fix(KM.DISPLACEMENT_Z)

        surface[6].Fix(KM.DISPLACEMENT_X)
        surface[6].Fix(KM.DISPLACEMENT_Y)
        surface[6].Fix(KM.DISPLACEMENT_Z)

        # clamp bending
        surface[0].Fix(KM.ROTATION_X)
        surface[0].Fix(KM.ROTATION_Y)
        surface[0].Fix(KM.ROTATION_Z)

        surface[3].Fix(KM.ROTATION_X)
        surface[3].Fix(KM.ROTATION_Y)
        surface[3].Fix(KM.ROTATION_Z)

        surface[6].Fix(KM.ROTATION_X)
        surface[6].Fix(KM.ROTATION_Y)
        surface[6].Fix(KM.ROTATION_Z)

        # apply neumann conditions
        prop = model_part.GetProperties()[2]

        model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [3], prop)
        surface[2].SetSolutionStepValue(SMA.POINT_LOAD_Z, -50)

        model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [6], prop)
        surface[5].SetSolutionStepValue(SMA.POINT_LOAD_Z, -100)

        model_part.CreateNewCondition('PointLoadCondition3D1N', 3, [9], prop)
        surface[8].SetSolutionStepValue(SMA.POINT_LOAD_Z, -50)

        # setup solver
        model_part.SetBufferSize(1)
        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        linear_solver = linear_solver_factory.ConstructSolver(
            KM.Parameters(r'{"solver_type": "LinearSolversApplication.sparse_lu"}'))

        relative_tolerance = 1e-8
        absolute_tolerance = 1e-7

        conv_criteria = KM.ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(0)

        maximum_iterations = 1
        compute_reactions = False
        reform_dofs_at_each_iteration = False
        move_mesh_flag = True

        solver = KM.ResidualBasedNewtonRaphsonStrategy(
            model_part,
            time_scheme,
            linear_solver,
            conv_criteria,
            maximum_iterations,
            compute_reactions,
            reform_dofs_at_each_iteration,
            move_mesh_flag
        )

        solver.SetEchoLevel(0)
        model_part.CloneTimeStep(1)
        solver.Solve()

        return surface

    def testCantileverOneQuadraticSpanWithoutParameterDistortion(self):
        def create_geometry(model_part):
            node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0)
            node2 = model_part.CreateNewNode(2, 2.5, 0.0, 0)
            node3 = model_part.CreateNewNode(3, 5.0, 0.0, 0)

            node4 = model_part.CreateNewNode(4, 0.0, 0.5, 0)
            node5 = model_part.CreateNewNode(5, 2.5, 0.5, 0)
            node6 = model_part.CreateNewNode(6, 5.0, 0.5, 0)

            node7 = model_part.CreateNewNode(7, 0.0, 1.0, 0)
            node8 = model_part.CreateNewNode(8, 2.5, 1.0, 0)
            node9 = model_part.CreateNewNode(9, 5.0, 1.0, 0)

            nodes = KM.NodesVector()
            nodes.append(node1)
            nodes.append(node2)
            nodes.append(node3)
            nodes.append(node4)
            nodes.append(node5)
            nodes.append(node6)
            nodes.append(node7)
            nodes.append(node8)
            nodes.append(node9)

            knots_u = KM.Vector(4)
            knots_u[0] = 0.0
            knots_u[1] = 0.0
            knots_u[2] = 5.0
            knots_u[3] = 5.0

            knots_v = KM.Vector(4)
            knots_v[0] = 0.0
            knots_v[1] = 0.0
            knots_v[2] = 1.0
            knots_v[3] = 1.0

            surface = KM.NurbsSurfaceGeometry3D(
                nodes, 2, 2, knots_u, knots_v)

            return surface

        surface = Shell6pElementTests.solve_cantilever(create_geometry)

        for node in surface:
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_X), 0)
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_Y), 0)

        self.assertAlmostEqual(surface[0].Z, 0.0)
        self.assertAlmostEqual(surface[3].Z, 0.0)
        self.assertAlmostEqual(surface[6].Z, 0.0)

        self.assertAlmostEqual(surface[1].Z, -0.000394798427724)
        self.assertAlmostEqual(surface[4].Z, -0.000286103781494)
        self.assertAlmostEqual(surface[7].Z, -0.000394798427727)

        self.assertAlmostEqual(surface[2].Z, -0.375659527565902)
        self.assertAlmostEqual(surface[5].Z, -0.375832346880606)
        self.assertAlmostEqual(surface[8].Z, -0.375659527565938)

    def testCantileverOneQuadraticSpanWithParameterDistortion(self):
        def create_geometry(model_part):
            node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0)
            node2 = model_part.CreateNewNode(2, 2.5, 0.0, 0)
            node3 = model_part.CreateNewNode(3, 5.0, 0.0, 0)

            node4 = model_part.CreateNewNode(4, 0.0, 0.5, 0)
            node5 = model_part.CreateNewNode(5, 2.5, 0.5, 0)
            node6 = model_part.CreateNewNode(6, 5.0, 0.5, 0)

            node7 = model_part.CreateNewNode(7, 0.0, 1.0, 0)
            node8 = model_part.CreateNewNode(8, 2.5, 1.0, 0)
            node9 = model_part.CreateNewNode(9, 5.0, 1.0, 0)

            nodes = KM.NodesVector()
            nodes.append(node1)
            nodes.append(node2)
            nodes.append(node3)
            nodes.append(node4)
            nodes.append(node5)
            nodes.append(node6)
            nodes.append(node7)
            nodes.append(node8)
            nodes.append(node9)

            knots_u = KM.Vector(4)
            knots_u[0] = 0.0
            knots_u[1] = 0.0
            knots_u[2] = 1.0
            knots_u[3] = 1.0

            knots_v = KM.Vector(4)
            knots_v[0] = 0.0
            knots_v[1] = 0.0
            knots_v[2] = 2.0
            knots_v[3] = 2.0

            surface = KM.NurbsSurfaceGeometry3D(
                nodes, 2, 2, knots_u, knots_v)

            return surface

        surface = Shell6pElementTests.solve_cantilever(create_geometry)

        for node in surface:
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_X), 0)
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_Y), 0)

        self.assertAlmostEqual(surface[0].Z, 0.0)
        self.assertAlmostEqual(surface[3].Z, 0.0)
        self.assertAlmostEqual(surface[6].Z, 0.0)

        self.assertAlmostEqual(surface[1].Z, -0.000394798427721)
        self.assertAlmostEqual(surface[4].Z, -0.000286103781492)
        self.assertAlmostEqual(surface[7].Z, -0.000394798427728)

        self.assertAlmostEqual(surface[2].Z, -0.375659527565157)
        self.assertAlmostEqual(surface[5].Z, -0.375832346879848)
        self.assertAlmostEqual(surface[8].Z, -0.375659527565172)


_GV = KM.KratosGlobals.GetVariable


def make_shell_model_part(model, name="IgaModelPart"):
    mp = model.CreateModelPart(name)
    mp.SetBufferSize(2)
    for var in (KM.DISPLACEMENT, KM.REACTION, KM.ROTATION, KM.REACTION_MOMENT,
                IGA.POINT_LOAD, IGA.DEAD_LOAD, SMA.POINT_LOAD):
        mp.AddNodalSolutionStepVariable(var)
    return mp


def make_surface_from_data(model_part, p_u, p_v, knots_u, knots_v, control_points, weights=None):
    """Build a (rational) NurbsSurfaceGeometry3D from hardcoded data. `control_points` in
    V-outer (u-inner) order; `knots_u/v` are the REDUCED Kratos knot vectors."""
    n_u = len(knots_u) - p_u + 1
    n_v = len(knots_v) - p_v + 1
    nodes = KM.NodesVector()
    nid = max((n.Id for n in model_part.Nodes), default=0) + 1
    for (x, y, z) in control_points:
        nodes.append(model_part.CreateNewNode(nid, float(x), float(y), float(z)))
        nid += 1
    ku = KM.Vector(len(knots_u))
    for i, v in enumerate(knots_u):
        ku[i] = float(v)
    kv = KM.Vector(len(knots_v))
    for i, v in enumerate(knots_v):
        kv[i] = float(v)
    if weights is not None:
        w = KM.Vector(len(weights))
        for i, v in enumerate(weights):
            w[i] = float(v)
        return KM.NurbsSurfaceGeometry3D(nodes, p_u, p_v, ku, kv, w), n_u, n_v
    return KM.NurbsSurfaceGeometry3D(nodes, p_u, p_v, ku, kv), n_u, n_v


def create_shell_elements(model_part, surface, element_name, props, deriv_order=3, eid_start=1):
    qpg = KM.GeometriesVector()
    surface.CreateQuadraturePointGeometries(qpg, deriv_order)
    for i in range(len(qpg)):
        model_part.CreateNewElement(element_name, eid_start + i, qpg[i], props)
    return qpg


def attach_isotropic_cl(props, E, nu, thickness, density=1.0):
    """Single-ply isotropic material (the shells build their own isotropic matrix from E, nu;
    the CL object sets the strain size)."""
    props.SetValue(KM.THICKNESS, thickness)
    props.SetValue(KM.YOUNG_MODULUS, E)
    props.SetValue(KM.POISSON_RATIO, nu)
    props.SetValue(KM.DENSITY, density)
    props.SetValue(KM.CONSTITUTIVE_LAW, SMA.LinearElasticPlaneStress2DLaw())


def apply_dead_load(main_mp, qp_geoms, force_vec, props_id=3, submp="DeadLoad"):
    """Uniform per-unit-area dead load (Scordelis self-weight) on every QP geometry."""
    sub = main_mp.CreateSubModelPart(submp) if not main_mp.HasSubModelPart(submp) else main_mp.GetSubModelPart(submp)
    if not main_mp.HasProperties(props_id):
        main_mp.CreateNewProperties(props_id)
    props = main_mp.GetProperties()[props_id]
    f = KM.Vector(3); f[0], f[1], f[2] = (float(v) for v in force_vec)
    cid = max((c.Id for c in main_mp.Conditions), default=0) + 1
    for i in range(len(qp_geoms)):
        cond = main_mp.CreateNewCondition("LoadCondition", cid + i, qp_geoms[i], props)
        cond.SetValue(IGA.DEAD_LOAD, f)
        sub.AddCondition(cond)


def add_point_load(main_mp, node, force_vec, submp="PointLoads", props_id=2):
    """PointLoadCondition3D1N at `node` carrying IGA.POINT_LOAD = force_vec."""
    sub = main_mp.GetSubModelPart(submp) if main_mp.HasSubModelPart(submp) else main_mp.CreateSubModelPart(submp)
    if not main_mp.HasProperties(props_id):
        main_mp.CreateNewProperties(props_id)
    props = main_mp.GetProperties()[props_id]
    cid = max((c.Id for c in main_mp.Conditions), default=0) + 1
    cond = main_mp.CreateNewCondition("PointLoadCondition3D1N", cid, [node.Id], props)
    f = KM.Vector(3); f[0], f[1], f[2] = (float(v) for v in force_vec)
    cond.SetValue(IGA.POINT_LOAD, f)
    sub.AddCondition(cond)
    if not sub.HasNode(node.Id):
        sub.AddNode(node, 0)
    return cond


_EDGE_GEOM_KEEPALIVE = []
_EDGE_ENDPOINTS = {
    "v=0": ((0.0, 0.0), (1.0, 0.0)),
    "v=1": ((0.0, 1.0), (1.0, 1.0)),
    "u=0": ((0.0, 0.0), (0.0, 1.0)),
    "u=1": ((1.0, 0.0), (1.0, 1.0)),
}


def _make_edge_curve_on_surface(surface, edge):
    p0, p1 = _EDGE_ENDPOINTS[edge]
    curve_nodes = KM.NodesVector()
    curve_nodes.append(KM.Node(1, p0[0], p0[1], 0.0))
    curve_nodes.append(KM.Node(2, p1[0], p1[1], 0.0))
    knots = KM.Vector(2); knots[0] = 0.0; knots[1] = 1.0
    curve_2d = KM.NurbsCurveGeometry2D(curve_nodes, 1, knots)
    return KM.NurbsCurveOnSurfaceGeometry3D(surface, curve_2d)


def add_rotation_penalty_on_edge(main_mp, surface, edge, penalty_factor,
                                 submp="RotPenalty", props_id=4, integration_order=3):
    """Weak zero-surface-normal-rotation BC along an isoparametric edge via
    SupportPenaltyRotationCondition on each edge integration point (Kirchhoff-Love)."""
    edge_curve = _make_edge_curve_on_surface(surface, edge)
    edge_qp = KM.GeometriesVector()
    edge_curve.CreateQuadraturePointGeometries(edge_qp, integration_order)
    sub = main_mp.GetSubModelPart(submp) if main_mp.HasSubModelPart(submp) else main_mp.CreateSubModelPart(submp)
    if not main_mp.HasProperties(props_id):
        main_mp.CreateNewProperties(props_id)
    props = main_mp.GetProperties()[props_id]
    props.SetValue(IGA.PENALTY_FACTOR, float(penalty_factor))
    cid = max((c.Id for c in main_mp.Conditions), default=0) + 1
    for i in range(len(edge_qp)):
        cond = main_mp.CreateNewCondition("SupportPenaltyRotationCondition", cid + i, edge_qp[i], props)
        sub.AddCondition(cond)
    _EDGE_GEOM_KEEPALIVE.append(edge_curve)
    _EDGE_GEOM_KEEPALIVE.append(edge_qp)
    return edge_qp


def add_dofs(main_mp, add_rotation_dofs=True):
    pairs = [(KM.DISPLACEMENT_X, KM.REACTION_X), (KM.DISPLACEMENT_Y, KM.REACTION_Y),
             (KM.DISPLACEMENT_Z, KM.REACTION_Z)]
    if add_rotation_dofs:
        pairs += [(KM.ROTATION_X, KM.REACTION_MOMENT_X), (KM.ROTATION_Y, KM.REACTION_MOMENT_Y),
                  (KM.ROTATION_Z, KM.REACTION_MOMENT_Z)]
    for d, r in pairs:
        KM.VariableUtils().AddDof(d, r, main_mp)


def solve_linear(main_mp, add_rotation_dofs=True):
    add_dofs(main_mp, add_rotation_dofs)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    linear_solver = KM.SkylineLUFactorizationSolver()
    bs = KM.ResidualBasedBlockBuilderAndSolver(linear_solver)
    bs.SetCalculateReactionsFlag(False)
    strat = KM.ResidualBasedLinearStrategy(main_mp, scheme, bs, False, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    main_mp.CloneTimeStep(1.0)
    strat.Solve()
    strat.FinalizeSolutionStep()
    return bool(strat.IsConverged())



#obstacle-course geometries
SCORDELIS = {
    "P_U": 3, "P_V": 3, "N_U": 7, "N_V": 7,
    "KNOTS_U": [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
    "KNOTS_V": [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
    "E": 432000000.0, "NU": 0.0, "T": 0.25, "RHO_G": 90.0, "DRILLING_PENALTY": 1e-05, "ROTATION_PENALTY": 100000000.0, "U_REF": -0.3024,
    "WEIGHTS": [
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.9899487701, 0.9899487701, 0.9899487701, 0.9899487701, 0.9899487701, 0.9899487701, 0.9899487701,
        0.9748719253, 0.9748719253, 0.9748719253, 0.9748719253, 0.9748719253, 0.9748719253, 0.9748719253,
        0.9673335029, 0.9673335029, 0.9673335029, 0.9673335029, 0.9673335029, 0.9673335029, 0.9673335029,
        0.9748719253, 0.9748719253, 0.9748719253, 0.9748719253, 0.9748719253, 0.9748719253, 0.9748719253,
        0.9899487701, 0.9899487701, 0.9899487701, 0.9899487701, 0.9899487701, 0.9899487701, 0.9899487701,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    ],
    "CPS": [
        [0.0, 0.0, 25.0],
        [0.0, 2.0833333333, 25.0],
        [0.0, 6.25, 25.0],
        [0.0, 12.5, 25.0],
        [0.0, 18.75, 25.0],
        [0.0, 22.9166666667, 25.0],
        [0.0, 25.0, 25.0],
        [1.4395532108, 0.0, 25.0],
        [1.4395532108, 2.0833333333, 25.0],
        [1.4395532108, 6.25, 25.0],
        [1.4395532108, 12.5, 25.0],
        [1.4395532108, 18.75, 25.0],
        [1.4395532108, 22.9166666667, 25.0],
        [1.4395532108, 25.0, 25.0],
        [4.3413705361, 0.0, 24.7500146442],
        [4.3413705361, 2.0833333333, 24.7500146442],
        [4.3413705361, 6.25, 24.7500146442],
        [4.3413705361, 12.5, 24.7500146442],
        [4.3413705361, 18.75, 24.7500146442],
        [4.3413705361, 22.9166666667, 24.7500146442],
        [4.3413705361, 25.0, 24.7500146442],
        [8.5949262547, 0.0, 23.6143658067],
        [8.5949262547, 2.0833333333, 23.6143658067],
        [8.5949262547, 6.25, 23.6143658067],
        [8.5949262547, 12.5, 23.6143658067],
        [8.5949262547, 18.75, 23.6143658067],
        [8.5949262547, 22.9166666667, 23.6143658067],
        [8.5949262547, 25.0, 23.6143658067],
        [12.5833199782, 0.0, 21.750190375],
        [12.5833199782, 2.0833333333, 21.750190375],
        [12.5833199782, 6.25, 21.750190375],
        [12.5833199782, 12.5, 21.750190375],
        [12.5833199782, 18.75, 21.750190375],
        [12.5833199782, 22.9166666667, 21.750190375],
        [12.5833199782, 25.0, 21.750190375],
        [14.9669285045, 0.0, 20.0764380453],
        [14.9669285045, 2.0833333333, 20.0764380453],
        [14.9669285045, 6.25, 20.0764380453],
        [14.9669285045, 12.5, 20.0764380453],
        [14.9669285045, 18.75, 20.0764380453],
        [14.9669285045, 22.9166666667, 20.0764380453],
        [14.9669285045, 25.0, 20.0764380453],
        [16.0696902422, 0.0, 19.151111078],
        [16.0696902422, 2.0833333333, 19.151111078],
        [16.0696902422, 6.25, 19.151111078],
        [16.0696902422, 12.5, 19.151111078],
        [16.0696902422, 18.75, 19.151111078],
        [16.0696902422, 22.9166666667, 19.151111078],
        [16.0696902422, 25.0, 19.151111078],
    ],
}

PINCHED = {
    "P_U": 5, "P_V": 5, "N_U": 11, "N_V": 11,
    "KNOTS_U": [0.0, 0.0, 0.0, 0.0, 0.0, 0.1666666667, 0.3333333333, 0.5, 0.6666666667, 0.8333333333, 1.0, 1.0, 1.0, 1.0, 1.0],
    "KNOTS_V": [0.0, 0.0, 0.0, 0.0, 0.0, 0.1666666667, 0.3333333333, 0.5, 0.6666666667, 0.8333333333, 1.0, 1.0, 1.0, 1.0, 1.0],
    "E": 3000000.0, "NU": 0.3, "T": 3.0, "P": 1.0, "DRILLING_PENALTY": 1e-05, "ROTATION_PENALTY": 100000000.0, "U_REF": -1.8248e-05,
    "WEIGHTS": [
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854,
        0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253,
        0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425,
        0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133,
        0.8454174678, 0.8454174678, 0.8454174678, 0.8454174678, 0.8454174678, 0.8454174678, 0.8454174678, 0.8454174678, 0.8454174678, 0.8454174678, 0.8454174678,
        0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133, 0.8616893133,
        0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425, 0.9007417425,
        0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253, 0.9446757253,
        0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854, 0.9804737854,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    ],
    "CPS": [
        [0.0, 0.0, 300.0],
        [0.0, 10.0, 300.0],
        [0.0, 30.0, 300.0],
        [0.0, 60.0, 300.0],
        [0.0, 100.0, 300.0],
        [0.0, 150.0, 300.0],
        [0.0, 200.0, 300.0],
        [0.0, 240.0, 300.0],
        [0.0, 270.0, 300.0],
        [0.0, 290.0, 300.0],
        [0.0, 300.0, 300.0],
        [14.4237773963, 0.0, 300.0],
        [14.4237773963, 10.0, 300.0],
        [14.4237773963, 30.0, 300.0],
        [14.4237773963, 60.0, 300.0],
        [14.4237773963, 100.0, 300.0],
        [14.4237773963, 150.0, 300.0],
        [14.4237773963, 200.0, 300.0],
        [14.4237773963, 240.0, 300.0],
        [14.4237773963, 270.0, 300.0],
        [14.4237773963, 290.0, 300.0],
        [14.4237773963, 300.0, 300.0],
        [44.1802936337, 0.0, 298.235726163],
        [44.1802936337, 10.0, 298.235726163],
        [44.1802936337, 30.0, 298.235726163],
        [44.1802936337, 60.0, 298.235726163],
        [44.1802936337, 100.0, 298.235726163],
        [44.1802936337, 150.0, 298.235726163],
        [44.1802936337, 200.0, 298.235726163],
        [44.1802936337, 240.0, 298.235726163],
        [44.1802936337, 270.0, 298.235726163],
        [44.1802936337, 290.0, 298.235726163],
        [44.1802936337, 300.0, 298.235726163],
        [89.9878980418, 0.0, 289.8232021078],
        [89.9878980418, 10.0, 289.8232021078],
        [89.9878980418, 30.0, 289.8232021078],
        [89.9878980418, 60.0, 289.8232021078],
        [89.9878980418, 100.0, 289.8232021078],
        [89.9878980418, 150.0, 289.8232021078],
        [89.9878980418, 200.0, 289.8232021078],
        [89.9878980418, 240.0, 289.8232021078],
        [89.9878980418, 270.0, 289.8232021078],
        [89.9878980418, 290.0, 289.8232021078],
        [89.9878980418, 300.0, 289.8232021078],
        [150.1006515141, 0.0, 266.1517600193],
        [150.1006515141, 10.0, 266.1517600193],
        [150.1006515141, 30.0, 266.1517600193],
        [150.1006515141, 60.0, 266.1517600193],
        [150.1006515141, 100.0, 266.1517600193],
        [150.1006515141, 150.0, 266.1517600193],
        [150.1006515141, 200.0, 266.1517600193],
        [150.1006515141, 240.0, 266.1517600193],
        [150.1006515141, 270.0, 266.1517600193],
        [150.1006515141, 290.0, 266.1517600193],
        [150.1006515141, 300.0, 266.1517600193],
        [216.2149635805, 0.0, 216.2149635805],
        [216.2149635805, 10.0, 216.2149635805],
        [216.2149635805, 30.0, 216.2149635805],
        [216.2149635805, 60.0, 216.2149635805],
        [216.2149635805, 100.0, 216.2149635805],
        [216.2149635805, 150.0, 216.2149635805],
        [216.2149635805, 200.0, 216.2149635805],
        [216.2149635805, 240.0, 216.2149635805],
        [216.2149635805, 270.0, 216.2149635805],
        [216.2149635805, 290.0, 216.2149635805],
        [216.2149635805, 300.0, 216.2149635805],
        [266.1517600193, 0.0, 150.1006515141],
        [266.1517600193, 10.0, 150.1006515141],
        [266.1517600193, 30.0, 150.1006515141],
        [266.1517600193, 60.0, 150.1006515141],
        [266.1517600193, 100.0, 150.1006515141],
        [266.1517600193, 150.0, 150.1006515141],
        [266.1517600193, 200.0, 150.1006515141],
        [266.1517600193, 240.0, 150.1006515141],
        [266.1517600193, 270.0, 150.1006515141],
        [266.1517600193, 290.0, 150.1006515141],
        [266.1517600193, 300.0, 150.1006515141],
        [289.8232021078, 0.0, 89.9878980418],
        [289.8232021078, 10.0, 89.9878980418],
        [289.8232021078, 30.0, 89.9878980418],
        [289.8232021078, 60.0, 89.9878980418],
        [289.8232021078, 100.0, 89.9878980418],
        [289.8232021078, 150.0, 89.9878980418],
        [289.8232021078, 200.0, 89.9878980418],
        [289.8232021078, 240.0, 89.9878980418],
        [289.8232021078, 270.0, 89.9878980418],
        [289.8232021078, 290.0, 89.9878980418],
        [289.8232021078, 300.0, 89.9878980418],
        [298.235726163, 0.0, 44.1802936337],
        [298.235726163, 10.0, 44.1802936337],
        [298.235726163, 30.0, 44.1802936337],
        [298.235726163, 60.0, 44.1802936337],
        [298.235726163, 100.0, 44.1802936337],
        [298.235726163, 150.0, 44.1802936337],
        [298.235726163, 200.0, 44.1802936337],
        [298.235726163, 240.0, 44.1802936337],
        [298.235726163, 270.0, 44.1802936337],
        [298.235726163, 290.0, 44.1802936337],
        [298.235726163, 300.0, 44.1802936337],
        [300.0, 0.0, 14.4237773963],
        [300.0, 10.0, 14.4237773963],
        [300.0, 30.0, 14.4237773963],
        [300.0, 60.0, 14.4237773963],
        [300.0, 100.0, 14.4237773963],
        [300.0, 150.0, 14.4237773963],
        [300.0, 200.0, 14.4237773963],
        [300.0, 240.0, 14.4237773963],
        [300.0, 270.0, 14.4237773963],
        [300.0, 290.0, 14.4237773963],
        [300.0, 300.0, 14.4237773963],
        [300.0, 0.0, 0.0],
        [300.0, 10.0, 0.0],
        [300.0, 30.0, 0.0],
        [300.0, 60.0, 0.0],
        [300.0, 100.0, 0.0],
        [300.0, 150.0, 0.0],
        [300.0, 200.0, 0.0],
        [300.0, 240.0, 0.0],
        [300.0, 270.0, 0.0],
        [300.0, 290.0, 0.0],
        [300.0, 300.0, 0.0],
    ],
}

HEMI = {
    "P_U": 4, "P_V": 4, "N_U": 10, "N_V": 10,
    "KNOTS_U": [0.0, 0.0, 0.0, 0.0, 0.1666666667, 0.3333333333, 0.5, 0.6666666667, 0.8333333333, 1.0, 1.0, 1.0, 1.0],
    "KNOTS_V": [0.0, 0.0, 0.0, 0.0, 0.1666666667, 0.3333333333, 0.5, 0.6666666667, 0.8333333333, 1.0, 1.0, 1.0, 1.0],
    "E": 68250000.0, "NU": 0.3, "T": 0.04, "LOAD": 1.0, "DRILLING_PENALTY": 1e-05, "ROTATION_PENALTY": 100000000.0, "U_REF": 0.0924,
    "WEIGHTS": [
        1.0, 0.9755922318, 0.9322006438, 0.8833851073, 0.8508414163, 0.8508414163, 0.8833851073, 0.9322006438, 0.9755922318, 1.0,
        0.9755922318, 0.9517802027, 0.9094477065, 0.8618236484, 0.8300742763, 0.8300742763, 0.8618236484, 0.9094477065, 0.9517802027, 0.9755922318,
        0.9322006438, 0.9094477065, 0.8689980403, 0.8234921658, 0.7931549161, 0.7931549161, 0.8234921658, 0.8689980403, 0.9094477065, 0.9322006438,
        0.8833851073, 0.8618236484, 0.8234921658, 0.7803692478, 0.7516206359, 0.7516206359, 0.7803692478, 0.8234921658, 0.8618236484, 0.8833851073,
        0.8508414163, 0.8300742763, 0.7931549161, 0.7516206359, 0.7239311158, 0.7239311158, 0.7516206359, 0.7931549161, 0.8300742763, 0.8508414163,
        0.8508414163, 0.8300742763, 0.7931549161, 0.7516206359, 0.7239311158, 0.7239311158, 0.7516206359, 0.7931549161, 0.8300742763, 0.8508414163,
        0.8833851073, 0.8618236484, 0.8234921658, 0.7803692478, 0.7516206359, 0.7516206359, 0.7803692478, 0.8234921658, 0.8618236484, 0.8833851073,
        0.9322006438, 0.9094477065, 0.8689980403, 0.8234921658, 0.7931549161, 0.7931549161, 0.8234921658, 0.8689980403, 0.9094477065, 0.9322006438,
        0.9755922318, 0.9517802027, 0.9094477065, 0.8618236484, 0.8300742763, 0.8300742763, 0.8618236484, 0.9094477065, 0.9517802027, 0.9755922318,
        1.0, 0.9755922318, 0.9322006438, 0.8833851073, 0.8508414163, 0.8508414163, 0.8833851073, 0.9322006438, 0.9755922318, 1.0,
    ],
    "CPS": [
        [0.0, 0.0, 10.0],
        [0.6039978915, 0.0, 10.0],
        [1.8551948627, 0.0, 9.9006731081],
        [3.7634682614, 0.0, 9.4235138729],
        [6.1367218686, 0.0, 8.0955671183],
        [8.0955671183, 0.0, 6.1367218686],
        [9.4235138729, 0.0, 3.7634682614],
        [9.9006731081, 0.0, 1.8551948627],
        [10.0, 0.0, 0.6039978915],
        [10.0, 0.0, 0.0],
        [0.0, 0.0, 10.0],
        [0.6039978915, 0.0364813453, 10.0],
        [1.8551948627, 0.1120533785, 9.9006731081],
        [3.7634682614, 0.2273126895, 9.4235138729],
        [6.1367218686, 0.370656707, 8.0955671183],
        [8.0955671183, 0.488970547, 6.1367218686],
        [9.4235138729, 0.569178251, 3.7634682614],
        [9.9006731081, 0.5979985682, 1.8551948627],
        [10.0, 0.6039978915, 0.6039978915],
        [10.0, 0.6039978915, 0.0],
        [0.0, 0.0, 10.0],
        [0.5979985682, 0.1120533785, 10.0],
        [1.8367677887, 0.3441747978, 9.9006731081],
        [3.7260869008, 0.6981966984, 9.4235138729],
        [6.0757677176, 1.1384814884, 8.0955671183],
        [8.0151563663, 1.5018854528, 6.1367218686],
        [9.3299130385, 1.7482454525, 3.7634682614],
        [9.8023327993, 1.8367677887, 1.8551948627],
        [9.9006731081, 1.8551948627, 0.6039978915],
        [9.9006731081, 1.8551948627, 0.0],
        [0.0, 0.0, 10.0],
        [0.569178251, 0.2273126895, 10.0],
        [1.7482454525, 0.6981966984, 9.9006731081],
        [3.5465095371, 1.4163693354, 9.4235138729],
        [5.7829483663, 2.3095357981, 8.0955671183],
        [7.6288689049, 3.0467409908, 6.1367218686],
        [8.8802613713, 3.5465095371, 3.7634682614],
        [9.3299130385, 3.7260869008, 1.8551948627],
        [9.4235138729, 3.7634682614, 0.6039978915],
        [9.4235138729, 3.7634682614, 0.0],
        [0.0, 0.0, 10.0],
        [0.488970547, 0.370656707, 10.0],
        [1.5018854528, 1.1384814884, 9.9006731081],
        [3.0467409908, 2.3095357981, 9.4235138729],
        [4.9680243774, 3.7659355293, 8.0955671183],
        [6.5538206967, 4.9680243774, 6.1367218686],
        [7.6288689049, 5.7829483663, 3.7634682614],
        [8.0151563663, 6.0757677176, 1.8551948627],
        [8.0955671183, 6.1367218686, 0.6039978915],
        [8.0955671183, 6.1367218686, 0.0],
        [0.0, 0.0, 10.0],
        [0.370656707, 0.488970547, 10.0],
        [1.1384814884, 1.5018854528, 9.9006731081],
        [2.3095357981, 3.0467409908, 9.4235138729],
        [3.7659355293, 4.9680243774, 8.0955671183],
        [4.9680243774, 6.5538206967, 6.1367218686],
        [5.7829483663, 7.6288689049, 3.7634682614],
        [6.0757677176, 8.0151563663, 1.8551948627],
        [6.1367218686, 8.0955671183, 0.6039978915],
        [6.1367218686, 8.0955671183, 0.0],
        [0.0, 0.0, 10.0],
        [0.2273126895, 0.569178251, 10.0],
        [0.6981966984, 1.7482454525, 9.9006731081],
        [1.4163693354, 3.5465095371, 9.4235138729],
        [2.3095357981, 5.7829483663, 8.0955671183],
        [3.0467409908, 7.6288689049, 6.1367218686],
        [3.5465095371, 8.8802613713, 3.7634682614],
        [3.7260869008, 9.3299130385, 1.8551948627],
        [3.7634682614, 9.4235138729, 0.6039978915],
        [3.7634682614, 9.4235138729, 0.0],
        [0.0, 0.0, 10.0],
        [0.1120533785, 0.5979985682, 10.0],
        [0.3441747978, 1.8367677887, 9.9006731081],
        [0.6981966984, 3.7260869008, 9.4235138729],
        [1.1384814884, 6.0757677176, 8.0955671183],
        [1.5018854528, 8.0151563663, 6.1367218686],
        [1.7482454525, 9.3299130385, 3.7634682614],
        [1.8367677887, 9.8023327993, 1.8551948627],
        [1.8551948627, 9.9006731081, 0.6039978915],
        [1.8551948627, 9.9006731081, 0.0],
        [0.0, 0.0, 10.0],
        [0.0364813453, 0.6039978915, 10.0],
        [0.1120533785, 1.8551948627, 9.9006731081],
        [0.2273126895, 3.7634682614, 9.4235138729],
        [0.370656707, 6.1367218686, 8.0955671183],
        [0.488970547, 8.0955671183, 6.1367218686],
        [0.569178251, 9.4235138729, 3.7634682614],
        [0.5979985682, 9.9006731081, 1.8551948627],
        [0.6039978915, 10.0, 0.6039978915],
        [0.6039978915, 10.0, 0.0],
        [0.0, 0.0, 10.0],
        [0.0, 0.6039978915, 10.0],
        [0.0, 1.8551948627, 9.9006731081],
        [0.0, 3.7634682614, 9.4235138729],
        [0.0, 6.1367218686, 8.0955671183],
        [0.0, 8.0955671183, 6.1367218686],
        [0.0, 9.4235138729, 3.7634682614],
        [0.0, 9.9006731081, 1.8551948627],
        [0.0, 10.0, 0.6039978915],
        [0.0, 10.0, 0.0],
    ],
}


#shell obstacle course for all three shells
def _is_kl(name):
    return name == "Shell3pElement"


def _obstacle_solve(g, element_name, disp_bcs, rot6_bcs, kl_edges, apply_load, monitor):
    model = KM.Model()
    mp = make_shell_model_part(model)
    surf, n_u, n_v = make_surface_from_data(
        mp, g["P_U"], g["P_V"], g["KNOTS_U"], g["KNOTS_V"], g["CPS"], g["WEIGHTS"])
    props = mp.GetProperties()[1]
    attach_isotropic_cl(props, g["E"], g["NU"], g["T"], density=g.get("RHO_G", 1.0))
    kl = _is_kl(element_name)
    if not kl:
        props.SetValue(_GV("DRILLING_PENALTY"), g["DRILLING_PENALTY"])
    qpg = create_shell_elements(mp, surf, element_name, props, deriv_order=4 if kl else 3)
    disp_bcs(surf, n_u, n_v)
    if kl:
        for edge in kl_edges:
            add_rotation_penalty_on_edge(mp, surf, edge, g["ROTATION_PENALTY"])
    else:
        rot6_bcs(surf, n_u, n_v)
    apply_load(g, mp, surf, n_u, n_v, qpg)
    solve_linear(mp, add_rotation_dofs=not kl)
    return monitor(g, surf, n_u, n_v)


def _scordelis_disp_bcs(surf, n_u, n_v):
    nd = lambda uu, vv: surf[vv * n_u + uu]
    for vv in range(n_v):
        nd(0, vv).Fix(KM.DISPLACEMENT_X); nd(0, vv).Fix(KM.DISPLACEMENT_Z)
    for vv in range(n_v): 
        nd(n_u - 1, vv).Fix(KM.DISPLACEMENT_Y)
    for uu in range(n_u):
        nd(uu, 0).Fix(KM.DISPLACEMENT_X)


def _scordelis_rot6_bcs(surf, n_u, n_v):
    nd = lambda uu, vv: surf[vv * n_u + uu]
    for vv in range(n_v):
        nd(0, vv).Fix(KM.ROTATION_Y)
    for vv in range(n_v):                                  
        nd(n_u - 1, vv).Fix(KM.ROTATION_X); nd(n_u - 1, vv).Fix(KM.ROTATION_Z)
    for uu in range(n_u):                                  
        nd(uu, 0).Fix(KM.ROTATION_Y); nd(uu, 0).Fix(KM.ROTATION_Z)


def _scordelis_load(g, mp, surf, n_u, n_v, qpg):
    apply_dead_load(mp, qpg, [0.0, 0.0, -g["RHO_G"]])

 # Punkt D, vertical u_z
def _scordelis_monitor(g, surf, n_u, n_v):
    return float(surf[(n_v - 1) * n_u + (n_u - 1)].GetSolutionStepValue(KM.DISPLACEMENT)[2])


def _solve_scordelis(name):
    return _obstacle_solve(SCORDELIS, name, _scordelis_disp_bcs, _scordelis_rot6_bcs,
                           ("u=1", "v=0"), _scordelis_load, _scordelis_monitor)


def _pinched_disp_bcs(surf, n_u, n_v):
    nd = lambda uu, vv: surf[vv * n_u + uu]
    for vv in range(n_v):
        nd(0, vv).Fix(KM.DISPLACEMENT_Y)
        nd(n_u - 1, vv).Fix(KM.DISPLACEMENT_X); nd(n_u - 1, vv).Fix(KM.DISPLACEMENT_Z)
    for uu in range(n_u):
        nd(uu, 0).Fix(KM.DISPLACEMENT_X)
        nd(uu, n_v - 1).Fix(KM.DISPLACEMENT_Z)


def _pinched_rot6_bcs(surf, n_u, n_v):
    nd = lambda uu, vv: surf[vv * n_u + uu]
    for vv in range(n_v):
        nd(0, vv).Fix(KM.ROTATION_X); nd(0, vv).Fix(KM.ROTATION_Z)
        nd(n_u - 1, vv).Fix(KM.ROTATION_Y)
    for uu in range(n_u):
        nd(uu, 0).Fix(KM.ROTATION_Y); nd(uu, 0).Fix(KM.ROTATION_Z)
        nd(uu, n_v - 1).Fix(KM.ROTATION_X); nd(uu, n_v - 1).Fix(KM.ROTATION_Y)


def _pinched_load(g, mp, surf, n_u, n_v, qpg):
    add_point_load(mp, surf[0], [0.0, 0.0, -g["P"] / 4.0])


def _pinched_monitor(g, surf, n_u, n_v):
    return float(surf[0].GetSolutionStepValue(KM.DISPLACEMENT)[2])


def _solve_pinched(name):
    return _obstacle_solve(PINCHED, name, _pinched_disp_bcs, _pinched_rot6_bcs,
                           ("v=0", "u=0", "v=1"), _pinched_load, _pinched_monitor)


def _hemi_disp_bcs(surf, n_u, n_v):
    nd = lambda uu, vv: surf[vv * n_u + uu]
    for uu in range(n_u):
        nd(uu, 0).Fix(KM.DISPLACEMENT_Y)
        nd(uu, n_v - 1).Fix(KM.DISPLACEMENT_X)


def _hemi_rot6_bcs(surf, n_u, n_v):
    nd = lambda uu, vv: surf[vv * n_u + uu]
    for uu in range(n_u):
        nd(uu, 0).Fix(KM.ROTATION_X); nd(uu, 0).Fix(KM.ROTATION_Z)
        nd(uu, n_v - 1).Fix(KM.ROTATION_Y); nd(uu, n_v - 1).Fix(KM.ROTATION_Z)


def _hemi_load(g, mp, surf, n_u, n_v, qpg):
    add_point_load(mp, surf[n_u - 1], [+g["LOAD"], 0.0, 0.0])
    add_point_load(mp, surf[(n_v - 1) * n_u + (n_u - 1)], [0.0, -g["LOAD"], 0.0])


def _hemi_monitor(g, surf, n_u, n_v):
    node_A = surf[n_u - 1]
    d = node_A.GetSolutionStepValue(KM.DISPLACEMENT)
    X, Y, Z = node_A.X, node_A.Y, node_A.Z
    rn = (X * X + Y * Y + Z * Z) ** 0.5 or 1.0
    return (d[0] * X + d[1] * Y + d[2] * Z) / rn


def _solve_hemisphere(name):
    return _obstacle_solve(HEMI, name, _hemi_disp_bcs, _hemi_rot6_bcs,
                           ("v=0", "v=1"), _hemi_load, _hemi_monitor)


_OBSTACLE = {
    "scordelis":  (_solve_scordelis,  SCORDELIS["U_REF"], 0.03),
    "pinched":    (_solve_pinched,    PINCHED["U_REF"],   0.05),
    "hemisphere": (_solve_hemisphere, HEMI["U_REF"],      0.05),
}


class IsotropicObstacleCourseShellTests(KratosUnittest.TestCase):
    """Scordelis-Lo roof, pinched cylinder and Belytschko hemisphere for the
    Shell3p (KL), Shell6p and Shell6p_Sommerwerk (RM) elements."""

    def _check(self, problem, element_name):
        solve, ref, rtol = _OBSTACLE[problem]
        val = solve(element_name)
        self.assertAlmostEqual(val / ref, 1.0, delta=rtol,
            msg=f"{element_name} {problem}: value {val:.6e} vs ref {ref:.6e} (ratio {val/ref:.4f})")

    def test_Shell3p_scordelis(self):            self._check("scordelis", "Shell3pElement")
    def test_Shell6p_scordelis(self):            self._check("scordelis", "Shell6pElement")
    def test_Shell6pSommerwerk_scordelis(self):  self._check("scordelis", "Shell6pElement_Sommerwerk")

    def test_Shell3p_pinched(self):              self._check("pinched", "Shell3pElement")
    def test_Shell6p_pinched(self):              self._check("pinched", "Shell6pElement")
    def test_Shell6pSommerwerk_pinched(self):    self._check("pinched", "Shell6pElement_Sommerwerk")

    def test_Shell3p_hemisphere(self):           self._check("hemisphere", "Shell3pElement")
    @unittest.expectedFailure
    def test_Shell6p_hemisphere(self):
        # KNOWN FAILURE: the legacy Shell6pElement locks on the Belytschko hemisphere
        # (ratio ~0.03); needs ANS/DSG. Marked expectedFailure so a future fix shows up
        # as an 'unexpected success'.
        self._check("hemisphere", "Shell6pElement")
    def test_Shell6pSommerwerk_hemisphere(self): self._check("hemisphere", "Shell6pElement_Sommerwerk")


if __name__ == "__main__":
    KratosUnittest.main()
