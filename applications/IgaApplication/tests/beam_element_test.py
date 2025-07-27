import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest
import numpy as np


class TestBeam4pElement(KratosUnittest.TestCase):

    def analytical_solution_bending(x, F, E, I, L):
        return F / (6 * E * I) * (3 * L * x**2 - x**3)
    
    def analytical_solution_torsion(x,M,E, It, nu = 0):
        G = 1 / (2 * (1 + nu)) * E
        return M / (G * It) * x
    
    def analytical_solution_normalforce(x,F,E,A):
        return (F / (E * A)) * x

    def solve_linear_examples(last_node_force, last_node_moment):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

        # create property for truss elements
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
        beam_properties.SetValue(IGA.CENTER_LINE_ROTATION, np.array([[0, 0], [0.25 , 0], [0.5 , 0],[0.75 , 0], [1, 0]]))
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

        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(quadrature_point_geometries,3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('IsogeometricBeamElement', element_id, quadrature_point_geometries[i], beam_properties)
            element_id += 1

        #Compute T0 and N0 vectors automatically using parent curve
        IGA.ComputeBeamVectorsProcess(model_part, curve).ExecuteInitialize()

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)

        # apply dirichlet conditions
        nodes[0].Fix(KM.DISPLACEMENT_X)
        nodes[0].Fix(KM.DISPLACEMENT_Y)
        nodes[0].Fix(KM.DISPLACEMENT_Z)
        nodes[0].Fix(KM.ROTATION_X)
        nodes[1].Fix(KM.DISPLACEMENT_Y)
        nodes[1].Fix(KM.DISPLACEMENT_Z)

        # apply neumann conditions
        prop = model_part.GetProperties()[2]
        force_condition= model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [nodes[len(controllpoints)-1].Id], prop)
        force_condition.SetValue(SMA.POINT_LOAD, last_node_force)
        moment_condition = model_part.CreateNewCondition('PointMomentCondition3D1N', 2, [nodes[len(controllpoints)-1].Id], prop)
        moment_condition.SetValue(SMA.POINT_MOMENT, last_node_moment )

        # setup solver
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
        beam_properties.SetValue(IGA.CENTER_LINE_ROTATION, np.array([[0, 0], [0.25 , 0], [0.5 , 0],[0.75 , 0], [1, 0]]))
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

        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(quadrature_point_geometries,3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('IsogeometricBeamElement', element_id, quadrature_point_geometries[i], beam_properties)
            element_id += 1

        #Compute T0 and N0 vectors automatically using parent curve
        IGA.ComputeBeamVectorsProcess(model_part, curve).ExecuteInitialize()

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)

        # apply dirichlet conditions
        nodes[0].Fix(KM.DISPLACEMENT_X)
        nodes[0].Fix(KM.DISPLACEMENT_Y)
        nodes[0].Fix(KM.DISPLACEMENT_Z)
        nodes[0].Fix(KM.ROTATION_X)
        nodes[1].Fix(KM.DISPLACEMENT_X)
        nodes[1].Fix(KM.DISPLACEMENT_Y)
        nodes[1].Fix(KM.DISPLACEMENT_Z)

        prop = model_part.GetProperties()[2]
        #force_condition= model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [nodes[len(controllpoints)-1].Id], prop)
        #force_condition.SetValue(SMA.POINT_LOAD,  KM.Array3([0,  0, 1]))

        force_condition= model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [nodes[2].Id], prop)
        force_condition.SetValue(SMA.POINT_LOAD,  KM.Array3([0,  0, load]))

        # setup solver
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
        beam_properties.SetValue(IGA.CENTER_LINE_ROTATION, np.array([[0, 0], [0.25 , 0], [0.5 , 0],[0.75 , 0], [1, 0]]))
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

        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(quadrature_point_geometries,3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('IsogeometricBeamElement', element_id, quadrature_point_geometries[i], beam_properties)
            element_id += 1

        #Compute T0 and N0 vectors automatically using parent curve
        IGA.ComputeBeamVectorsProcess(model_part, curve).ExecuteInitialize()
        beam_properties.SetValue(IGA.T_0,  KM.Array3([-1,  0, 0])) #[-0.95281728851577396, 0.020230986128019116, 0.30286947998393016]
        beam_properties.SetValue(IGA.N_0, KM.Array3([0,  0, -1])) #[-0.30286947998393016, -0, -0.95281728851577396]

        
        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)

        # apply dirichlet conditions
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

        # setup solver
        model_part.SetBufferSize(1)

        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()

        linear_solver = linear_solver_factory.ConstructSolver(KM.Parameters(
            r'{"solver_type": "skyline_lu_factorization"}'))

        relative_tolerance = 1e-7
        absolute_tolerance = 1e-7

        conv_criteria = KM.ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(0)

        maximum_iterations = 20
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

    # def testClampedFzNonlinear(self):
    #     last_node= TestBeam4pElement.solve_nonlinear_example2()

    #     disp = np.array(last_node.GetSolutionStepValue(KM.DISPLACEMENT)) #[np.float64(-0.13288290406702538), np.float64(-0.29808295882441466), np.float64(-1.1774740489689273)]
    #     print(disp)
    #     pass 



