from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import new_linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class TrussDiscreteElementTest(KratosUnittest.TestCase):

    def solve(create_geometry):
        model_part = ModelPart('Model')

        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(REACTION)
        model_part.AddNodalSolutionStepVariable(POINT_LOAD)

        # create property for truss elements

        truss_properties = model_part.GetProperties()[1]
        truss_properties.SetValue(CROSS_AREA      , 0.01  )
        truss_properties.SetValue(PRESTRESS_CAUCHY, 0     )
        truss_properties.SetValue(YOUNG_MODULUS   , 210000)
        truss_properties.SetValue(POISSON_RATIO   , 0     )
        truss_properties.SetValue(DENSITY         , 7856  )

        # create a node based geometry

        node_1, node_2, curve = create_geometry(model_part)

        # create elements for each integration point

        integration_points = [integration_point for span in curve.Spans()
            for integration_point in IntegrationPoints.Points1D(curve.Degree + 1,
            span)]

        shapes = CurveShapeEvaluator(Degree=curve.Degree, Order=1)

        for i, (t, weight) in enumerate(integration_points):
            shapes.Compute(curve.Knots(), t)

            node_indices = [curve.Node(j).Id for j in
                range(shapes.FirstNonzeroPole, shapes.LastNonzeroPole + 1)]

            element = model_part.CreateNewElement('TrussDiscreteElement', i + 1,
                node_indices, truss_properties)

            n_0 = Vector(shapes.NumberOfNonzeroPoles)
            for i in range(shapes.NumberOfNonzeroPoles):
                n_0[i] = shapes(0, i)

            n_1 = Matrix(shapes.NumberOfNonzeroPoles, 1)
            for i in range(shapes.NumberOfNonzeroPoles):
                n_1[i, 0] = shapes(1, i)

            element.SetValue(INTEGRATION_WEIGHT, weight)
            element.SetValue(SHAPE_FUNCTION_VALUES, n_0)
            element.SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, n_1)

        # add dofs

        VariableUtils().AddDof(DISPLACEMENT_X, REACTION_X, model_part)
        VariableUtils().AddDof(DISPLACEMENT_Y, REACTION_Y, model_part)
        VariableUtils().AddDof(DISPLACEMENT_Z, REACTION_Z, model_part)

        # apply dirichlet conditions

        node_1.Fix(DISPLACEMENT_X)
        node_1.Fix(DISPLACEMENT_Y)
        node_1.Fix(DISPLACEMENT_Z)

        node_2.Fix(DISPLACEMENT_Y)
        node_2.Fix(DISPLACEMENT_Z)

        # apply neumann conditions

        prop = model_part.GetProperties()[2]

        model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [node_2.Id],
            prop)

        node_2.SetSolutionStepValue(POINT_LOAD_X, 1000)

        # setup solver

        model_part.SetBufferSize(1)

        time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
            r'{"solver_type": "SkylineLUFactorizationSolver"}'))

        relative_tolerance = 1e-7
        absolute_tolerance = 1e-7

        conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(0)

        maximum_iterations = 100
        compute_reactions = False
        reform_dofs_at_each_iteration = True
        move_mesh_flag = True

        solver = ResidualBasedNewtonRaphsonStrategy(
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

        return node_1, node_2, curve

    def testOneLinearSpanWithoutParameterDistortion(self):
        def create_geometry(model_part):
            node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            node_2 = model_part.CreateNewNode(2, 2.0, 0.0, 0.0)

            node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
            node_2.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)

            curve = NodeCurveGeometry3D(1, 2)

            curve.SetKnot(Index=0, Value=0.0)
            curve.SetKnot(Index=1, Value=2.0)

            curve.SetNode(Index=0, Value=node_1)
            curve.SetNode(Index=1, Value=node_2)

            return node_1, node_2, curve

        node_1, node_2, _ = TrussDiscreteElementTest.solve(create_geometry)

        self.assertAlmostEqual(node_1.X, 0.0               )
        self.assertAlmostEqual(node_1.Y, 0.0               )
        self.assertAlmostEqual(node_1.Z, 0.0               )

        self.assertAlmostEqual(node_2.X, 2.6268671884353934)
        self.assertAlmostEqual(node_2.Y, 0.0               )
        self.assertAlmostEqual(node_2.Z, 0.0               )

    def testOneLinearSpanWithParameterDistortion(self):
        def create_geometry(model_part):
            node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            node_2 = model_part.CreateNewNode(2, 2.0, 0.0, 0.0)

            node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
            node_2.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)

            curve = NodeCurveGeometry3D(1, 2)

            curve.SetKnot(Index=0, Value=0.0)
            curve.SetKnot(Index=1, Value=5.0)

            curve.SetNode(Index=0, Value=node_1)
            curve.SetNode(Index=1, Value=node_2)

            return node_1, node_2, curve

        node_1, node_2, _ = TrussDiscreteElementTest.solve(create_geometry)

        self.assertAlmostEqual(node_1.X, 0.0               )
        self.assertAlmostEqual(node_1.Y, 0.0               )
        self.assertAlmostEqual(node_1.Z, 0.0               )

        self.assertAlmostEqual(node_2.X, 2.6268671884353934)
        self.assertAlmostEqual(node_2.Y, 0.0               )
        self.assertAlmostEqual(node_2.Z, 0.0               )

    def testTwoLinearSpans(self):
        def create_geometry(model_part):
            node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
            node_3 = model_part.CreateNewNode(3, 2.0, 0.0, 0.0)

            node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
            node_2.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
            node_3.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)

            curve = NodeCurveGeometry3D(1, 3)

            curve.SetKnot(Index=0, Value=0.0)
            curve.SetKnot(Index=1, Value=1.0)
            curve.SetKnot(Index=2, Value=2.0)

            curve.SetNode(Index=0, Value=node_1)
            curve.SetNode(Index=1, Value=node_2)
            curve.SetNode(Index=2, Value=node_3)

            return node_1, node_3, curve

        node_1, node_2, _ = TrussDiscreteElementTest.solve(create_geometry)

        self.assertAlmostEqual(node_1.X, 0.0               )
        self.assertAlmostEqual(node_1.Y, 0.0               )
        self.assertAlmostEqual(node_1.Z, 0.0               )

        self.assertAlmostEqual(node_2.X, 2.6268671884353934)
        self.assertAlmostEqual(node_2.Y, 0.0               )
        self.assertAlmostEqual(node_2.Z, 0.0               )

    def testTwoLinearSpansWithParameterDistortion(self):
        def create_geometry(model_part):
            node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
            node_3 = model_part.CreateNewNode(3, 2.0, 0.0, 0.0)

            node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
            node_2.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
            node_3.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)

            curve = NodeCurveGeometry3D(1, 3)

            curve.SetKnot(Index=0, Value=0.0)
            curve.SetKnot(Index=1, Value=3.0)
            curve.SetKnot(Index=2, Value=4.0)

            curve.SetNode(Index=0, Value=node_1)
            curve.SetNode(Index=1, Value=node_2)
            curve.SetNode(Index=2, Value=node_3)

            return node_1, node_3, curve

        node_1, node_2, _ = TrussDiscreteElementTest.solve(create_geometry)

        self.assertAlmostEqual(node_1.X, 0.0               )
        self.assertAlmostEqual(node_1.Y, 0.0               )
        self.assertAlmostEqual(node_1.Z, 0.0               )

        self.assertAlmostEqual(node_2.X, 2.6268671884353934)
        self.assertAlmostEqual(node_2.Y, 0.0               )
        self.assertAlmostEqual(node_2.Z, 0.0               )

    def testOneQuadraticSpan(self):
        def create_geometry(model_part):
            node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
            node_3 = model_part.CreateNewNode(3, 2.0, 0.0, 0.0)

            node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
            node_2.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
            node_3.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)

            curve = NodeCurveGeometry3D(2, 3)

            curve.SetKnot(Index=0, Value=0.0)
            curve.SetKnot(Index=1, Value=0.0)
            curve.SetKnot(Index=2, Value=2.0)
            curve.SetKnot(Index=3, Value=2.0)

            curve.SetNode(Index=0, Value=node_1)
            curve.SetNode(Index=1, Value=node_2)
            curve.SetNode(Index=2, Value=node_3)

            return node_1, node_3, curve

        node_1, node_2, _ = TrussDiscreteElementTest.solve(create_geometry)

        self.assertAlmostEqual(node_1.X, 0.0               )
        self.assertAlmostEqual(node_1.Y, 0.0               )
        self.assertAlmostEqual(node_1.Z, 0.0               )

        self.assertAlmostEqual(node_2.X, 2.6268671884353934)
        self.assertAlmostEqual(node_2.Y, 0.0               )
        self.assertAlmostEqual(node_2.Z, 0.0               )
