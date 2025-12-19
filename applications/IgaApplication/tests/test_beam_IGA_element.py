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


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class BeamIGAElementTests(KratosUnittest.TestCase):

    def solve_cantilever_thin_beam_element(create_geometry):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

        # create property for beam elements

        beam_properties = model_part.GetProperties()[1]
        beam_properties.SetValue(KM.THICKNESS, 1.0)
        beam_properties.SetValue(SMA.CROSS_AREA, 1.0)
        beam_properties.SetValue(KM.YOUNG_MODULUS, 200000)
        beam_properties.SetValue(KM.POISSON_RATIO, 0)
        beam_properties.SetValue(KM.CONSTITUTIVE_LAW, SMA.BeamConstitutiveLaw())

        # create a nurbs curve
        curve = create_geometry(model_part)

        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()

        curve.CreateQuadraturePointGeometries(quadrature_point_geometries, 3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('BeamThinElement2D', element_id, quadrature_point_geometries[i], beam_properties)
            element_id += 1

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)

        # apply dirichlet conditions
        curve[0].Fix(KM.DISPLACEMENT_X)
        curve[0].Fix(KM.DISPLACEMENT_Y)

        # clamp bending
        curve[1].Fix(KM.DISPLACEMENT_Y)

        # apply neumann conditions
        prop = model_part.GetProperties()[2]

        model_part.CreateNewCondition('PointLoadCondition2D1N', 1, [4], prop)
        curve[3].SetSolutionStepValue(SMA.POINT_LOAD_Y, -1.0)

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

        return curve

    def solve_cantilever_thick_beam_element(create_geometry):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)
        model_part.AddNodalSolutionStepVariable(IGA.CROSS_SECTIONAL_ROTATION)

        # create property for beam elements

        beam_properties = model_part.GetProperties()[1]
        beam_properties.SetValue(KM.THICKNESS, 1.0)
        beam_properties.SetValue(SMA.CROSS_AREA, 1.0)
        beam_properties.SetValue(KM.YOUNG_MODULUS, 200000)
        beam_properties.SetValue(KM.POISSON_RATIO, 0)
        beam_properties.SetValue(KM.CONSTITUTIVE_LAW, SMA.BeamConstitutiveLaw())

        # create a nurbs curve
        curve = create_geometry(model_part)

        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()

        curve.CreateQuadraturePointGeometries(quadrature_point_geometries, 3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('BeamThickElement2D', element_id, quadrature_point_geometries[i], beam_properties)
            element_id += 1

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(IGA.CROSS_SECTIONAL_ROTATION, model_part)

        # apply dirichlet conditions
        curve[0].Fix(KM.DISPLACEMENT_X)
        curve[0].Fix(KM.DISPLACEMENT_Y)
        curve[0].Fix(IGA.CROSS_SECTIONAL_ROTATION)

        # apply neumann conditions
        prop = model_part.GetProperties()[2]

        model_part.CreateNewCondition('PointLoadCondition2D1N', 1, [4], prop)
        curve[3].SetSolutionStepValue(SMA.POINT_LOAD_Y, -1.0)

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

        solver.SetEchoLevel(3)
        model_part.CloneTimeStep(1)
        solver.Solve()

        return curve

    def testCantileverThinBeamTwoQuadraticSpan(self):
        def create_geometry(model_part):
            node1 = model_part.CreateNewNode(1, 0.0,  0.0, 0)
            node2 = model_part.CreateNewNode(2, 2.5, 0.0, 0)
            node3 = model_part.CreateNewNode(3, 7.5, 0.0, 0)
            node4 = model_part.CreateNewNode(4, 10.0,  0.0, 0)

            nodes = KM.NodesVector()
            nodes.append(node1)
            nodes.append(node2)
            nodes.append(node3)
            nodes.append(node4)

            knots_u = KM.Vector(5)
            knots_u[0] = 0.0
            knots_u[1] = 0.0
            knots_u[2] = 0.5
            knots_u[3] = 1.0
            knots_u[4] = 1.0

            curve = KM.NurbsCurveGeometry2D(nodes, 2, knots_u)

            return curve

        # Thin beam test
        curve_1 = BeamIGAElementTests.solve_cantilever_thin_beam_element(create_geometry)

        for node in curve_1:
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_X), 0)

        self.assertAlmostEqual(curve_1[0].Y, 0.0)
        self.assertAlmostEqual(curve_1[1].Y, 0.0) 
        self.assertAlmostEqual(curve_1[2].Y, -0.01125)
        self.assertAlmostEqual(curve_1[3].Y, -0.01875)

        # Thick beam test
        curve_2 = BeamIGAElementTests.solve_cantilever_thick_beam_element(create_geometry)

        for node in curve_2:
            KratosUnittest.TestCase().assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_X), 0)

        KratosUnittest.TestCase().assertAlmostEqual(curve_2[0].Y, 0.0)
        KratosUnittest.TestCase().assertAlmostEqual(curve_2[1].Y, -0.000235479)
        KratosUnittest.TestCase().assertAlmostEqual(curve_2[2].Y, -0.0106722)
        KratosUnittest.TestCase().assertAlmostEqual(curve_2[3].Y, -0.0184077)