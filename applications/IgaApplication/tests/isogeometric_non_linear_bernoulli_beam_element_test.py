import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest
import numpy as np

class TestBeam4pElement(KratosUnittest.TestCase):

    def analytical_solution_bending(x, F, E, I, L):
        return F / (6 * E * I) * (3 * L * x**2 - x**3)

    def analytical_solution_torsion(x, M, E, It, nu = 0):
        G = 1 / (2 * (1 + nu)) * E
        return M / (G * It) * x

    def analytical_solution_normalforce(x, F, E, A):
        return (F / (E * A)) * x

    def solve_linear_examples(last_node_force, last_node_moment):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

        w = 2.5
        h = 5
        beam_properties = model_part.GetProperties()[0]
        beam_properties.SetValue(IGA.CROSS_AREA, w*h)
        beam_properties.SetValue(IGA.HEIGHT, h)
        beam_properties.SetValue(IGA.WIDTH, w)
        beam_properties.SetValue(IGA.I_T,  w*h**3/3)
        beam_properties.SetValue(IGA.I_N, w * h**3 / 12)
        beam_properties.SetValue(IGA.I_V, h * w**3 / 12)
        beam_properties.SetValue(KM.YOUNG_MODULUS   , 70000)
        beam_properties.SetValue(KM.POISSON_RATIO   , 0)
        beam_properties.SetValue(KM.DENSITY         , 7856)
        beam_properties.SetValue(IGA.LOCAL_AXIS_ORIENTATION, np.array([[0, 0, 0, 1], [0.25, 0, 0, 1], [0.5, 0, 0, 1], [0.75, 0, 0, 1], [1, 0, 0, 1]]))
        beam_properties.SetValue(KM.CONSTITUTIVE_LAW,IGA.BernoulliBeamElasticConstitutiveLaw())
        controllpoints = [[0.0, 0.0, 0.0],[1.0, 0.0, 0.0],[2.0, 0.0, 0.0],[3.0, 0.0, 0.0],[4.0, 0.0, 0.0]]
        knotvector = [0,0,0,0,1,1,1,1]

        nodes = KM.NodesVector()
        for i, ctrlpt in enumerate(controllpoints, 1):
            node = model_part.CreateNewNode(i, *ctrlpt)
            nodes.append(node)

        knots = KM.Vector(len(knotvector))
        for i, knot in enumerate(knotvector):
            knots[i] = knot

        curve = KM.NurbsCurveGeometry3D(nodes, 4, knots)

        quadrature_point_geometries = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(quadrature_point_geometries,3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('IsogeometricBeamElement', element_id, quadrature_point_geometries[i], beam_properties)
            element_id += 1

        IGA.ComputeBeamVectorsProcess(model_part, curve).ExecuteInitialize()

        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)

        nodes[0].Fix(KM.DISPLACEMENT_X)
        nodes[0].Fix(KM.DISPLACEMENT_Y)
        nodes[0].Fix(KM.DISPLACEMENT_Z)
        nodes[0].Fix(KM.ROTATION_X)
        nodes[1].Fix(KM.DISPLACEMENT_Y)
        nodes[1].Fix(KM.DISPLACEMENT_Z)

        prop = model_part.GetProperties()[2]
        force_condition= model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [nodes[len(controllpoints)-1].Id], prop)
        force_condition.SetValue(SMA.POINT_LOAD, last_node_force)
        moment_condition = model_part.CreateNewCondition('PointMomentCondition3D1N', 2, [nodes[len(controllpoints)-1].Id], prop)
        moment_condition.SetValue(SMA.POINT_MOMENT, last_node_moment )

        model_part.SetBufferSize(1)

        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()

        linear_solver = linear_solver_factory.ConstructSolver(KM.Parameters(
            r'{"solver_type": "LinearSolversApplication.sparse_lu"}'))

        relative_tolerance = 1e-7
        absolute_tolerance = 1e-7

        conv_criteria = KM.ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(0)

        maximum_iterations = 1
        compute_reactions = False
        reform_dofs_at_each_iteration = False
        move_mesh_flag = False

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

        return nodes, curve

    def solve_nonlinear_example1():
        load = -15

        controllpoints = [[0.0, 0.0, 0.0], [-0.60412500000000002, 0.0, 0.19230000000000003], [-1.1560208333333331, -0.047874999999999987, 0.36798749999999997], [-2.1115562499999996, -0.19984583333333331, 0.67512760416666673], [-2.9186499999999995, -0.41159583333333322, 0.94091666666666662], [-3.5380499999999997, -0.62918489583333348, 1.1495924479166666], [-4.139524999999999, -0.88037447916666678, 1.3563434895833333], [-4.7245499999999998, -1.1637291666666667, 1.5627458333333333], [-5.287725, -1.4909666666666668, 1.7699208333333334], [-5.8283447916666677, -1.8611583333333335, 1.9781567708333334], [-6.3503989583333347, -2.2607291666666667, 2.1869088541666666], [-6.849591666666667, -2.6874000000000002, 2.3957541666666664], [-7.313341666666668, -3.1501874999999999, 2.6043916666666669], [-7.740093229166666, -3.6477776041666661, 2.8127442708333334], [-8.1393161458333321, -4.1685255208333327, 3.0209588541666665], [-8.5090833333333329, -4.7094541666666672, 3.2291499999999997], [-8.8360708333333342, -5.2762541666666669, 3.4373999999999993], [-9.1194135416666668, -5.8673156249999998, 3.6457296874999994], [-9.370705208333332, -6.4737281250000009, 3.8540984374999989], [-9.5884125000000004, -7.0920291666666673, 4.0624749999999992], [-9.7578750000000021, -7.7242041666666674, 4.2708374999999998], [-9.8785942708333359, -8.3684583333333329, 4.4791796874999994], [-9.9642338541666682, -9.0192166666666651, 4.687510937499999], [-10.014333333333333, -9.672762500000001, 4.8958375000000007], [-10.014308333333332, -10.327237500000003, 5.1041625000000002], [-9.9641640625000001, -10.980783333333333, 5.3124880208333334], [-9.8784953124999983, -11.631541666666667, 5.5208151041666671], [-9.7577833333333324, -12.275795833333333, 5.7291458333333329], [-9.5883958333333332, -12.907970833333334, 5.9374833333333346], [-9.3708125000000013, -13.526271874999999, 6.1458270833333337], [-9.1196250000000028, -14.132684375, 6.3541729166666672], [-8.8362708333333337, -14.723745833333334, 6.5625166666666663], [-8.5090333333333348, -15.290545833333333, 6.7708541666666662], [-8.1388416666666679, -15.831474479166666, 6.9791848958333329], [-7.7392708333333342, -16.352222395833333, 7.1875119791666666], [-7.3126000000000007, -16.849812499999999, 7.3958375000000007], [-6.8498125000000014, -17.3126, 7.6041624999999993], [-6.3522223958333344, -17.739270833333332, 7.8124880208333334], [-5.8314744791666673, -18.138841666666668, 8.0208151041666653], [-5.2905458333333346, -18.509033333333335, 8.2291458333333338], [-4.723745833333334, -18.83627083333333, 8.4374833333333346], [-4.1326843750000011, -19.119625520833328, 8.6458270833333337], [-3.3746687500000014, -19.433612499999995, 8.9062593749999976], [-2.3873527777777794, -19.752688888888891, 9.2361347222222214], [-1.3089500000000003, -19.942550000000001, 9.5833500000000011], [0.0, -20.057449999999999, 10.0], [1.3089499999999998, -19.942550000000001, 10.416649999999999], [2.6030333333333333, -19.714716666666664, 10.833308333333331], [3.83765, -19.264949999999999, 11.249999999999996], [5.0286833333333343, -18.709958333333333, 11.666691666666667], [6.1048499999999999, -17.956, 12.083350000000003], [7.1116750000000009, -17.111674999999998, 12.5], [7.9559999999999995, -16.104850000000003, 12.916649999999999], [8.7099583333333346, -15.028683333333333, 13.333300000000001], [9.2649499999999989, -13.837650000000002, 13.74995], [9.7147083333333342, -12.603033333333332, 14.166608333333336], [9.942499999999999, -11.308949999999999, 14.583300000000001], [10.057383333333334, -10.0, 14.999991666666666], [9.9426000000000005, -8.6910500000000006, 15.416649999999999], [9.7148250000000012, -7.3969666666666676, 15.833324999999999], [9.2647499999999976, -6.16235, 16.250099999999996], [8.7096583333333335, -4.9713166666666657, 16.666799999999999], [7.9567999999999994, -3.8951499999999992, 17.083100000000002], [7.1127666666666656, -2.8883249999999996, 17.499650000000003], [6.1017999999999999, -2.044, 17.917599999999997], [5.0245499999999996, -1.2900416666666668, 18.334583333333327], [3.8491000000000009, -0.73504999999999987, 18.746299999999998], [2.6186499999999997, -0.28528333333333333, 19.161662499999998], [1.2664, -0.057450000000000001, 19.596875000000001], [0.60412500000000002, 0.0, 19.807700000000001], [0.0, 0.0, 20.0]]
        knotvector = [0.0, 0.0, 0.0, 0.041666666666666664, 0.083333333333333287, 0.10416666666666664, 0.125, 0.14583333333333331, 0.16666666666666666, 0.1875, 0.20833333333333331, 0.22916666666666666, 0.25, 0.27083333333333331, 0.29166666666666663, 0.3125, 0.33333333333333331, 0.35416666666666663, 0.375, 0.39583333333333331, 0.41666666666666663, 0.4375, 0.45833333333333331, 0.47916666666666663, 0.5, 0.52083333333333337, 0.54166666666666674, 0.5625, 0.58333333333333337, 0.60416666666666674, 0.625, 0.64583333333333337, 0.66666666666666674, 0.6875, 0.70833333333333337, 0.72916666666666674, 0.75, 0.77083333333333337, 0.79166666666666674, 0.8125, 0.83333333333333337, 0.85416666666666674, 0.875, 0.89583333333333326, 0.91666666666666663, 0.95833333333333337, 1.0, 1.0, 1.0, 13.0, 14.0, 14.0, 15.0, 15.0, 16.0, 16.0, 17.0, 17.0, 18.0, 18.0, 19.0, 19.0, 20.0, 20.0, 21.0, 21.0, 22.0, 22.0, 23.0, 23.0, 24.0, 24.0, 24.0, 24.0]

        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

        w = 0.1
        h = 0.1
        beam_properties = model_part.GetProperties()[0]
        beam_properties.SetValue(IGA.CROSS_AREA, w*h)
        beam_properties.SetValue(IGA.HEIGHT, h)
        beam_properties.SetValue(IGA.WIDTH, w)
        beam_properties.SetValue(IGA.I_T,  w*h**3/3)
        beam_properties.SetValue(IGA.I_V, w * h**3 / 12)
        beam_properties.SetValue(IGA.I_N, h * w**3 / 12)
        beam_properties.SetValue(KM.YOUNG_MODULUS   , 1e+10)
        beam_properties.SetValue(KM.POISSON_RATIO   , 0)
        beam_properties.SetValue(KM.DENSITY         , 1)
        beam_properties.SetValue(IGA.LOCAL_AXIS_ORIENTATION, np.array([[0, 0, 0, 1], [0.25, 0, 0, 1], [0.5, 0, 0, 1], [0.75, 0, 0, 1], [1, 0, 0, 1]]))
        beam_properties.SetValue(KM.CONSTITUTIVE_LAW,IGA.BernoulliBeamElasticConstitutiveLaw())

        nodes = KM.NodesVector()
        for i, ctrlpt in enumerate(controllpoints, 1):
            node = model_part.CreateNewNode(i, *ctrlpt)
            nodes.append(node)
        last_node = node
        knots = KM.Vector(len(knotvector))
        for i, knot in enumerate(knotvector):
            knots[i] = knot

        curve = KM.NurbsCurveGeometry3D(nodes, 4, knots)

        quadrature_point_geometries = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(quadrature_point_geometries,3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('IsogeometricBeamElement', element_id, quadrature_point_geometries[i], beam_properties)
            element_id += 1

        IGA.ComputeBeamVectorsProcess(model_part, curve).ExecuteInitialize()

        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)

        nodes[0].Fix(KM.DISPLACEMENT_X)
        nodes[0].Fix(KM.DISPLACEMENT_Y)
        nodes[0].Fix(KM.DISPLACEMENT_Z)
        nodes[0].Fix(KM.ROTATION_X)
        nodes[1].Fix(KM.DISPLACEMENT_X)
        nodes[1].Fix(KM.DISPLACEMENT_Y)
        nodes[1].Fix(KM.DISPLACEMENT_Z)

        prop = model_part.GetProperties()[2]

        force_condition= model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [nodes[2].Id], prop)
        force_condition.SetValue(SMA.POINT_LOAD,  KM.Array3([0,  0, load]))

        model_part.SetBufferSize(1)
        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        linear_solver = linear_solver_factory.ConstructSolver(
            KM.Parameters(r'{"solver_type": "LinearSolversApplication.sparse_lu"}'))

        relative_tolerance = 1e-8
        absolute_tolerance = 1e-7

        conv_criteria = KM.ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(0)

        maximum_iterations = 10
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

        solver.SetEchoLevel(1)
        model_part.CloneTimeStep(1)
        solver.Solve()

        return last_node

    def solve_nonlinear_example2():
        k  = 0.6
        load = k * 25000
        controllpoints = [[0.0, 0.0, 0.0], [-0.60412500000000002, 0.0, 0.19230000000000003], [-1.1560208333333331, -0.047874999999999987, 0.36798749999999997], [-2.1115562499999996, -0.19984583333333331, 0.67512760416666673], [-2.9186499999999995, -0.41159583333333322, 0.94091666666666662], [-3.5380499999999997, -0.62918489583333348, 1.1495924479166666], [-4.139524999999999, -0.88037447916666678, 1.3563434895833333], [-4.7245499999999998, -1.1637291666666667, 1.5627458333333333], [-5.287725, -1.4909666666666668, 1.7699208333333334], [-5.8283447916666677, -1.8611583333333335, 1.9781567708333334], [-6.3503989583333347, -2.2607291666666667, 2.1869088541666666], [-6.849591666666667, -2.6874000000000002, 2.3957541666666664], [-7.313341666666668, -3.1501874999999999, 2.6043916666666669], [-7.740093229166666, -3.6477776041666661, 2.8127442708333334], [-8.1393161458333321, -4.1685255208333327, 3.0209588541666665], [-8.5090833333333329, -4.7094541666666672, 3.2291499999999997], [-8.8360708333333342, -5.2762541666666669, 3.4373999999999993], [-9.1194135416666668, -5.8673156249999998, 3.6457296874999994], [-9.370705208333332, -6.4737281250000009, 3.8540984374999989], [-9.5884125000000004, -7.0920291666666673, 4.0624749999999992], [-9.7578750000000021, -7.7242041666666674, 4.2708374999999998], [-9.8785942708333359, -8.3684583333333329, 4.4791796874999994], [-9.9642338541666682, -9.0192166666666651, 4.687510937499999], [-10.014333333333333, -9.672762500000001, 4.8958375000000007], [-10.014308333333332, -10.327237500000003, 5.1041625000000002], [-9.9641640625000001, -10.980783333333333, 5.3124880208333334], [-9.8784953124999983, -11.631541666666667, 5.5208151041666671], [-9.7577833333333324, -12.275795833333333, 5.7291458333333329], [-9.5883958333333332, -12.907970833333334, 5.9374833333333346], [-9.3708125000000013, -13.526271874999999, 6.1458270833333337], [-9.1196250000000028, -14.132684375, 6.3541729166666672], [-8.8362708333333337, -14.723745833333334, 6.5625166666666663], [-8.5090333333333348, -15.290545833333333, 6.7708541666666662], [-8.1388416666666679, -15.831474479166666, 6.9791848958333329], [-7.7392708333333342, -16.352222395833333, 7.1875119791666666], [-7.3126000000000007, -16.849812499999999, 7.3958375000000007], [-6.8498125000000014, -17.3126, 7.6041624999999993], [-6.3522223958333344, -17.739270833333332, 7.8124880208333334], [-5.8314744791666673, -18.138841666666668, 8.0208151041666653], [-5.2905458333333346, -18.509033333333335, 8.2291458333333338], [-4.723745833333334, -18.83627083333333, 8.4374833333333346], [-4.1326843750000011, -19.119625520833328, 8.6458270833333337], [-3.3746687500000014, -19.433612499999995, 8.9062593749999976], [-2.3873527777777794, -19.752688888888891, 9.2361347222222214], [-1.3089500000000003, -19.942550000000001, 9.5833500000000011], [0.0, -20.057449999999999, 10.0], [1.3089499999999998, -19.942550000000001, 10.416649999999999], [2.6030333333333333, -19.714716666666664, 10.833308333333331], [3.83765, -19.264949999999999, 11.249999999999996], [5.0286833333333343, -18.709958333333333, 11.666691666666667], [6.1048499999999999, -17.956, 12.083350000000003], [7.1116750000000009, -17.111674999999998, 12.5], [7.9559999999999995, -16.104850000000003, 12.916649999999999], [8.7099583333333346, -15.028683333333333, 13.333300000000001], [9.2649499999999989, -13.837650000000002, 13.74995], [9.7147083333333342, -12.603033333333332, 14.166608333333336], [9.942499999999999, -11.308949999999999, 14.583300000000001], [10.057383333333334, -10.0, 14.999991666666666], [9.9426000000000005, -8.6910500000000006, 15.416649999999999], [9.7148250000000012, -7.3969666666666676, 15.833324999999999], [9.2647499999999976, -6.16235, 16.250099999999996], [8.7096583333333335, -4.9713166666666657, 16.666799999999999], [7.9567999999999994, -3.8951499999999992, 17.083100000000002], [7.1127666666666656, -2.8883249999999996, 17.499650000000003], [6.1017999999999999, -2.044, 17.917599999999997], [5.0245499999999996, -1.2900416666666668, 18.334583333333327], [3.8491000000000009, -0.73504999999999987, 18.746299999999998], [2.6186499999999997, -0.28528333333333333, 19.161662499999998], [1.2664, -0.057450000000000001, 19.596875000000001], [0.60412500000000002, 0.0, 19.807700000000001], [0.0, 0.0, 20.0]]
        knotvector = [0.0, 0.0, 0.0, 0.041666666666666664, 0.083333333333333287, 0.10416666666666664, 0.125, 0.14583333333333331, 0.16666666666666666, 0.1875, 0.20833333333333331, 0.22916666666666666, 0.25, 0.27083333333333331, 0.29166666666666663, 0.3125, 0.33333333333333331, 0.35416666666666663, 0.375, 0.39583333333333331, 0.41666666666666663, 0.4375, 0.45833333333333331, 0.47916666666666663, 0.5, 0.52083333333333337, 0.54166666666666674, 0.5625, 0.58333333333333337, 0.60416666666666674, 0.625, 0.64583333333333337, 0.66666666666666674, 0.6875, 0.70833333333333337, 0.72916666666666674, 0.75, 0.77083333333333337, 0.79166666666666674, 0.8125, 0.83333333333333337, 0.85416666666666674, 0.875, 0.89583333333333326, 0.91666666666666663, 0.95833333333333337, 1.0, 1.0, 1.0, 13.0, 14.0, 14.0, 15.0, 15.0, 16.0, 16.0, 17.0, 17.0, 18.0, 18.0, 19.0, 19.0, 20.0, 20.0, 21.0, 21.0, 22.0, 22.0, 23.0, 23.0, 24.0, 24.0, 24.0, 24.0]

        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

        r = 0.5*0.1
        beam_properties = model_part.GetProperties()[0]
        beam_properties.SetValue(IGA.CROSS_AREA, np.pi * r**2)
        beam_properties.SetValue(IGA.HEIGHT, 2*r)
        beam_properties.SetValue(IGA.WIDTH, 2*r)
        beam_properties.SetValue(IGA.I_T,  np.pi * r**4/2)
        beam_properties.SetValue(IGA.I_V, np.pi * r**4/4)
        beam_properties.SetValue(IGA.I_N, np.pi * r**4/4)
        beam_properties.SetValue(KM.YOUNG_MODULUS   , 1e+10)
        beam_properties.SetValue(KM.POISSON_RATIO   , 0)
        beam_properties.SetValue(KM.DENSITY         , 0)
        beam_properties.SetValue(IGA.LOCAL_AXIS_ORIENTATION, np.array([[0, 0, 0, 1], [0.25, 0, 0, 1], [0.5, 0, 0, 1], [0.75, 0, 0, 1], [1, 0, 0, 1]]))
        beam_properties.SetValue(KM.CONSTITUTIVE_LAW,IGA.BernoulliBeamElasticConstitutiveLaw())

        nodes = KM.NodesVector()
        for i, ctrlpt in enumerate(controllpoints, 1):
            node = model_part.CreateNewNode(i, *ctrlpt)
            nodes.append(node)
        last_node = node
        knots = KM.Vector(len(knotvector))
        for i, knot in enumerate(knotvector):
            knots[i] = knot

        curve = KM.NurbsCurveGeometry3D(nodes, 4, knots)

        quadrature_point_geometries = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(quadrature_point_geometries,3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('IsogeometricBeamElement', element_id, quadrature_point_geometries[i], beam_properties)
            element_id += 1

        IGA.ComputeBeamVectorsProcess(model_part, curve).ExecuteInitialize()

        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)

        nodes[0].Fix(KM.DISPLACEMENT_X)
        nodes[0].Fix(KM.DISPLACEMENT_Y)
        nodes[0].Fix(KM.DISPLACEMENT_Z)
        nodes[0].Fix(KM.ROTATION_X)
        nodes[1].Fix(KM.DISPLACEMENT_X)
        nodes[1].Fix(KM.DISPLACEMENT_Y)
        nodes[1].Fix(KM.DISPLACEMENT_Z)

        prop = model_part.GetProperties()[2]
        force_condition= model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [nodes[2].Id], prop)
        force_condition.SetValue(SMA.POINT_LOAD,  KM.Array3([0,  0, -load]))

        model_part.SetBufferSize(1)

        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()

        linear_solver = linear_solver_factory.ConstructSolver(KM.Parameters(
            r'{"solver_type": "skyline_lu_factorization"}'))

        relative_tolerance = 1e-7
        absolute_tolerance = 1e-7

        conv_criteria = KM.ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(0)

        maximum_iterations = 3
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

        solver.SetEchoLevel(1)
        model_part.CloneTimeStep(1)
        solver.Solve()

        return last_node

    def testClampedFXLinear(self):
        last_node_force = KM.Array3([1,  0, 0])
        last_node_moment = KM.Array3([0,  0, 0])
        nodes, _ = TestBeam4pElement.solve_linear_examples(last_node_force, last_node_moment)

        self.assertAlmostEqual(nodes[0].X, 0.0               )
        self.assertAlmostEqual(nodes[0].Y, 0.0               )
        self.assertAlmostEqual(nodes[0].Z, 0.0               )

        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.DISPLACEMENT))[0], TestBeam4pElement.analytical_solution_normalforce(nodes[4].X, 1, 70000, 2.5*5))
        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.DISPLACEMENT))[1], 0.0               )
        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.DISPLACEMENT))[2], 0.0               )

    def testClampedFYLinear(self):
        last_node_force = KM.Array3([0,  -1, 0])
        last_node_moment = KM.Array3([0,  0, 0])
        nodes , _ = TestBeam4pElement.solve_linear_examples(last_node_force, last_node_moment)

        self.assertAlmostEqual(nodes[0].X, 0.0               )
        self.assertAlmostEqual(nodes[0].Y, 0.0               )
        self.assertAlmostEqual(nodes[0].Z, 0.0               )

        self.assertAlmostEqual(nodes[1].X, 1.0               )
        self.assertAlmostEqual(nodes[1].Y, 0.0               )
        self.assertAlmostEqual(nodes[1].Z, 0.0               )

        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.DISPLACEMENT))[0], 0.0               )
        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.DISPLACEMENT))[1], TestBeam4pElement.analytical_solution_bending(nodes[4].X, -1, 70000,5**3*2.5/12,4))
        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.DISPLACEMENT))[2], 0.0)

    def testClampedFZLinear(self):
        last_node_force = KM.Array3([0,  0, -1])
        last_node_moment = KM.Array3([0,  0, 0])
        nodes , _ = TestBeam4pElement.solve_linear_examples(last_node_force, last_node_moment)

        self.assertAlmostEqual(nodes[0].X, 0.0               )
        self.assertAlmostEqual(nodes[0].Y, 0.0               )
        self.assertAlmostEqual(nodes[0].Z, 0.0               )

        self.assertAlmostEqual(nodes[1].X, 1.0               )
        self.assertAlmostEqual(nodes[1].Y, 0.0               )
        self.assertAlmostEqual(nodes[1].Z, 0.0               )

        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.DISPLACEMENT))[0], 0.0               )
        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.DISPLACEMENT))[1], 0.0               )
        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.DISPLACEMENT))[2], TestBeam4pElement.analytical_solution_bending(nodes[4].X,-1, 70000,2.5**3*5/12,4))

    def testClampedMXLinear(self):
        last_node_force = KM.Array3([0,  0, 0])
        last_node_moment = KM.Array3([1,  0, 0])
        nodes, _ = TestBeam4pElement.solve_linear_examples(last_node_force, last_node_moment)

        self.assertAlmostEqual(nodes[0].GetSolutionStepValue(KM.ROTATION)[0], 0.0               )
        self.assertAlmostEqual(nodes[0].GetSolutionStepValue(KM.ROTATION)[1], 0.0               )
        self.assertAlmostEqual(nodes[0].GetSolutionStepValue(KM.ROTATION)[2], 0.0               )

        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.ROTATION))[0], TestBeam4pElement.analytical_solution_torsion(nodes[4].X, 1, 70000,2.5*5**3/3))
        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.ROTATION))[1], 0.0               )
        self.assertAlmostEqual(np.array(nodes[4].GetSolutionStepValue(KM.ROTATION))[2], 0.0               )

E_MOD = 70000.0
NU = 0.0
WIDTH = 2.5
HEIGHT = 5.0
LENGTH = 4.0
I_N_ = WIDTH * HEIGHT ** 3 / 12.0
I_V_ = HEIGHT * WIDTH ** 3 / 12.0
I_T_ = WIDTH * HEIGHT ** 3 / 3.0
G_MOD = E_MOD / (2.0 * (1.0 + NU))
DOFS = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z, KM.ROTATION_X)

class TestBeam4pElementNonlinear(KratosUnittest.TestCase):

    @staticmethod
    def build_cantilever(n_cp=5, degree=4):

        model = KM.Model()
        mp = model.CreateModelPart('IgaModelPart')
        for var in (KM.DISPLACEMENT, KM.ROTATION, KM.REACTION,
                    KM.REACTION_MOMENT, SMA.POINT_LOAD):
            mp.AddNodalSolutionStepVariable(var)

        props = mp.GetProperties()[0]
        props.SetValue(IGA.CROSS_AREA, WIDTH * HEIGHT)
        props.SetValue(IGA.HEIGHT, HEIGHT)
        props.SetValue(IGA.WIDTH, WIDTH)
        props.SetValue(IGA.I_T, I_T_)
        props.SetValue(IGA.I_N, I_N_)
        props.SetValue(IGA.I_V, I_V_)
        props.SetValue(KM.YOUNG_MODULUS, E_MOD)
        props.SetValue(KM.POISSON_RATIO, NU)
        props.SetValue(KM.DENSITY, 7856)
        props.SetValue(IGA.LOCAL_AXIS_ORIENTATION,
                       np.array([[0, 0, 0, 1], [0.25, 0, 0, 1], [0.5, 0, 0, 1],
                                 [0.75, 0, 0, 1], [1, 0, 0, 1]]))
        props.SetValue(KM.CONSTITUTIVE_LAW, IGA.BernoulliBeamElasticConstitutiveLaw())

        nodes = KM.NodesVector()
        for i in range(n_cp):
            nodes.append(mp.CreateNewNode(i + 1, LENGTH * i / (n_cp - 1), 0.0, 0.0))

        kv = ([0.0] * degree
              + [float(i) / (n_cp - degree) for i in range(1, n_cp - degree)]
              + [1.0] * degree)
        knots = KM.Vector(len(kv))
        for i, k in enumerate(kv):
            knots[i] = k

        curve = KM.NurbsCurveGeometry3D(nodes, degree, knots)
        curve.SetId(2)
        mp.AddGeometry(curve)

        from KratosMultiphysics.modeler_factory import KratosModelerFactory
        modeler_parameters = KM.Parameters(r'''[{
            "modeler_name": "RefinementModeler",
            "Parameters": {
                "echo_level": 0,
                "refinements_file_name": "modeler_tests/curve_h_refinements.iga.json"
            }
        }]''')
        refinement_modeler = KratosModelerFactory().ConstructListOfModelers(
            model, modeler_parameters)[0]
        refinement_modeler.PrepareGeometryModel()

        curve = mp.GetGeometry(2)
        vtk_curve = KM.BrepCurve(curve)
        vtk_curve.SetId(3)
        mp.AddGeometry(vtk_curve)
        nodes = KM.NodesVector()
        for node in mp.Nodes:
            nodes.append(node)

        qpg = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(qpg, 3)
        for eid in range(len(qpg)):
            mp.CreateNewElement('IsogeometricBeamElement', eid + 1, qpg[eid], props)

        IGA.ComputeBeamVectorsProcess(mp, curve).ExecuteInitialize()

        for dof, reaction in zip(DOFS, (KM.REACTION_X, KM.REACTION_Y,
                                        KM.REACTION_Z, KM.REACTION_MOMENT_X)):
            KM.VariableUtils().AddDof(dof, reaction, mp)

        for var in DOFS:
            nodes[0].Fix(var)
        nodes[1].Fix(KM.DISPLACEMENT_Y)
        nodes[1].Fix(KM.DISPLACEMENT_Z)

        mp.SetBufferSize(1)
        mp.CloneTimeStep(1)
        for element in mp.Elements:
            element.Initialize(mp.ProcessInfo)
        return model, mp, nodes

    @staticmethod
    def free_dofs(nodes):
        fixed = {4 * i + j
                 for i in range(len(nodes))
                 for j, var in enumerate(DOFS) if nodes[i].IsFixed(var)}
        return [i for i in range(4 * len(nodes)) if i not in fixed]

    @staticmethod
    def element_dof_index(element):
        idx = []
        for node in element.GetGeometry():
            base = 4 * (node.Id - 1)
            idx += [base, base + 1, base + 2, base + 3]
        return np.array(idx)

    @classmethod
    def apply_state(cls, nodes, u):
        for i in range(len(nodes)):
            nodes[i].SetSolutionStepValue(
                KM.DISPLACEMENT, KM.Array3([u[4 * i], u[4 * i + 1], u[4 * i + 2]]))
            nodes[i].SetSolutionStepValue(KM.ROTATION_X, float(u[4 * i + 3]))


    @staticmethod
    def write_vtk_output(model, mp, output_file_name="beam_deformation"):
        from KratosMultiphysics.IgaApplication.iga_vtk_output_process import IgaVTKOutputProcess

        settings = KM.Parameters(r'''{
            "model_part_name": "IgaModelPart",
            "output_file_name": "beam_deformation",
            "brep_curve_ids": [3],
            "nodal_solution_step_data_variables": ["DISPLACEMENT", "ROTATION"],
            "output_refinement_curve": [10],
            "output_control_type": "step",
            "output_frequency": 1
        }''')
        settings["output_file_name"].SetString(output_file_name)

        output_process = IgaVTKOutputProcess(model, settings)
        output_process.ExecuteInitialize()
        output_process.ExecuteBeforeSolutionLoop()
        mp.ProcessInfo[KM.TIME] = 1.0
        if output_process.IsOutputStep():
            output_process.PrintOutput()

    @classmethod
    def assemble(cls, mp, nodes, u, build_level=0):

        cls.apply_state(nodes, u)
        mp.ProcessInfo[IGA.BUILD_LEVEL] = build_level
        ndof = 4 * len(nodes)
        K = np.zeros((ndof, ndof))
        R = np.zeros(ndof)
        for element in mp.Elements:
            lhs, rhs = KM.Matrix(), KM.Vector()
            element.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            idx = cls.element_dof_index(element)
            K[np.ix_(idx, idx)] += np.array(lhs, copy=True)
            R[idx] += np.array(rhs, copy=True)
        return K, R

    @classmethod
    def newton(cls, mp, nodes, f_ext, n_steps=8, tol=1e-10, max_it=40):

        free = cls.free_dofs(nodes)
        u = np.zeros(4 * len(nodes))
        norms = []
        for step in range(1, n_steps + 1):
            f = f_ext * step / n_steps
            scale = max(1.0, np.linalg.norm(f[free]))
            norms = []
            for _ in range(max_it):
                K, R = cls.assemble(mp, nodes, u)
                residual = R + f
                norm = np.linalg.norm(residual[free])
                norms.append(norm)
                if norm < tol * scale:
                    break
                u[free] += np.linalg.solve(K[np.ix_(free, free)], residual[free])
        return u, norms

    def testTangentMatchesFiniteDifferences(self):
        model, mp, nodes = TestBeam4pElementNonlinear.build_cantilever()
        n_cp = len(nodes)
        ndof = 4 * n_cp

        u0 = np.zeros(ndof)
        for i in range(n_cp):
            x = nodes[i].X0
            u0[4 * i] = 0.004 * x
            u0[4 * i + 1] = 0.010 * x ** 2
            u0[4 * i + 2] = -0.006 * x ** 2
            u0[4 * i + 3] = 0.003 * x

        K, _ = TestBeam4pElementNonlinear.assemble(mp, nodes, u0)

        h = 1e-7
        K_fd = np.zeros_like(K)
        for j in range(ndof):
            up = u0.copy(); up[j] += h
            um = u0.copy(); um[j] -= h
            _, rp = TestBeam4pElementNonlinear.assemble(mp, nodes, up)
            _, rm = TestBeam4pElementNonlinear.assemble(mp, nodes, um)
            K_fd[:, j] = -(rp - rm) / (2.0 * h)

        error = np.abs(K - K_fd).max() / np.abs(K).max()
        self.assertLess(error, 1e-5,
                        'tangent stiffness is not the derivative of the internal '
                        'force (relative error %.3e) -- a missing or wrong '
                        'geometric-stiffness term' % error)

    def testQuadraticNewtonConvergence(self):
        model, mp, nodes = TestBeam4pElementNonlinear.build_cantilever()
        n_cp = len(nodes)
        load = 0.5 * E_MOD * I_N_ / LENGTH ** 2
        f_ext = np.zeros(4 * n_cp)
        f_ext[4 * (n_cp - 1) + 1] = -load

        _, norms = TestBeam4pElementNonlinear.newton(mp, nodes, f_ext, n_steps=4)

        self.assertGreaterEqual(len(norms), 3)
        self.assertLess(norms[-1], 1e-6, 'Newton did not converge')

        self.assertLess(norms[-1], norms[-2] ** 1.5,
                        'convergence is not quadratic: %s'
                        % ' '.join('%.2e' % n for n in norms))

    def testLargeDeflectionElastica(self):
        EI = E_MOD * I_N_
        alpha = 0.5
        load = alpha * EI / LENGTH ** 2
        ux_exact, uy_exact = -0.0636758499, -0.6485743026

        model, mp, nodes = TestBeam4pElementNonlinear.build_cantilever()
        n_cp = len(nodes)
        f_ext = np.zeros(4 * n_cp)
        f_ext[4 * (n_cp - 1) + 1] = -load
        u, _ = TestBeam4pElementNonlinear.newton(mp, nodes, f_ext, n_steps=8)
        TestBeam4pElementNonlinear.write_vtk_output(model, mp)

        uy = u[4 * (n_cp - 1) + 1]
        ux = u[4 * (n_cp - 1)]

        self.assertLess(abs(uy), abs(-load * LENGTH ** 3 / (3 * EI)) * 1.02)
        self.assertLess(abs(uy - uy_exact) / abs(uy_exact), 0.06)
        self.assertLess(abs(ux - ux_exact) / abs(ux_exact), 0.15)

    @classmethod
    def buckling_loads(cls, mp, nodes, load_vector, p_ref):
        n_cp = len(nodes)
        f_ext = np.zeros(4 * n_cp)
        for j, value in enumerate(load_vector):
            f_ext[4 * (n_cp - 1) + j] = value

        u, _ = cls.newton(mp, nodes, f_ext, n_steps=1)

        cls.apply_state(nodes, u)
        k_mat, _ = cls.assemble(mp, nodes, u, build_level=1)
        k_geo, _ = cls.assemble(mp, nodes, u, build_level=2)

        free = cls.free_dofs(nodes)
        A = k_mat[np.ix_(free, free)]
        B = k_geo[np.ix_(free, free)]

        mu = np.linalg.eigvals(-np.linalg.solve(A, B))
        mu = mu[np.abs(mu.imag) < 1e-8 * (1.0 + np.abs(mu.real))].real
        lam = np.sort(1.0 / mu[np.abs(mu) > 1e-12])
        return np.sort(lam[lam > 1e-8]) * p_ref

    def testEulerBucklingWeakAxis(self):
        p_ref = 1000.0
        model, mp, nodes = TestBeam4pElementNonlinear.build_cantilever()
        loads = TestBeam4pElementNonlinear.buckling_loads(
            mp, nodes, [-p_ref, 0.0, 0.0], p_ref)
        euler = np.pi ** 2 * E_MOD * I_V_ / (4.0 * LENGTH ** 2)
        self.assertAlmostEqual(loads[0] / euler, 1.0, delta=0.005)

    def testEulerBucklingStrongAxis(self):
        p_ref = 1000.0
        model, mp, nodes = TestBeam4pElementNonlinear.build_cantilever()
        loads = TestBeam4pElementNonlinear.buckling_loads(
            mp, nodes, [-p_ref, 0.0, 0.0], p_ref)
        euler = np.pi ** 2 * E_MOD * I_N_ / (4.0 * LENGTH ** 2)
        self.assertAlmostEqual(loads[1] / euler, 1.0, delta=0.005)

    def testLateralTorsionalBuckling(self):
        p_ref = 1000.0
        model, mp, nodes = TestBeam4pElementNonlinear.build_cantilever()
        loads = TestBeam4pElementNonlinear.buckling_loads(
            mp, nodes, [0.0, -p_ref, 0.0], p_ref)
        timoshenko = 4.013 * np.sqrt(E_MOD * I_V_ * G_MOD * I_T_) / LENGTH ** 2

        self.assertAlmostEqual(loads[0] / timoshenko, 1.0, delta=0.05)

    @KratosUnittest.expectedFailure
    def testLinearBendingIsMeshIndependent(self):
        EI = E_MOD * I_N_
        load = 1e-4 * EI / LENGTH ** 2
        reference = -load * LENGTH ** 3 / (3.0 * EI)
        for n_cp in (5, 7, 9, 13):
            model, mp, nodes = TestBeam4pElementNonlinear.build_cantilever(n_cp)
            f_ext = np.zeros(4 * n_cp)
            f_ext[4 * (n_cp - 1) + 1] = -load
            u, _ = TestBeam4pElementNonlinear.newton(mp, nodes, f_ext, n_steps=1)
            uy = u[4 * (n_cp - 1) + 1]
            self.assertLess(abs(uy - reference) / abs(reference), 1e-6,
                            'n_cp=%d (%d knot spans) is off by %.3e'
                            % (n_cp, mp.NumberOfElements() // 5,
                               abs(uy - reference) / abs(reference)))

if __name__ == "__main__":
    KratosUnittest.main()
