import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class TrussElementTests(KratosUnittest.TestCase):

    def solve(create_geometry):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

        # create property for truss elements
        truss_properties = model_part.GetProperties()[1]
        truss_properties.SetValue(IGA.CROSS_AREA      , 0.01  )
        truss_properties.SetValue(IGA.TRUSS_PRESTRESS_CAUCHY, 0     )
        truss_properties.SetValue(KM.YOUNG_MODULUS   , 210000)
        truss_properties.SetValue(KM.POISSON_RATIO   , 0     )
        truss_properties.SetValue(KM.DENSITY         , 7856  )
        truss_properties.SetValue(KM.CONSTITUTIVE_LAW,
            SMA.TrussConstitutiveLaw())

        # create a node based geometry
        node_1, node_2, curve = create_geometry(model_part)

        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(quadrature_point_geometries, 2)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('TrussElement', element_id, quadrature_point_geometries[i], truss_properties)
            element_id += 1

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

        # apply dirichlet conditions
        node_1.Fix(KM.DISPLACEMENT_X)
        node_1.Fix(KM.DISPLACEMENT_Y)
        node_1.Fix(KM.DISPLACEMENT_Z)

        node_2.Fix(KM.DISPLACEMENT_Y)
        node_2.Fix(KM.DISPLACEMENT_Z)

        # apply neumann conditions
        prop = model_part.GetProperties()[2]
        model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [node_2.Id], prop)
        node_2.SetSolutionStepValue(SMA.POINT_LOAD_X, 1000)

        # setup solver
        model_part.SetBufferSize(1)

        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()

        linear_solver = linear_solver_factory.ConstructSolver(KM.Parameters(
            r'{"solver_type": "skyline_lu_factorization"}'))

        relative_tolerance = 1e-7
        absolute_tolerance = 1e-7

        conv_criteria = KM.ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(0)

        maximum_iterations = 100
        compute_reactions = True
        reform_dofs_at_each_iteration = True
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

        solver.SetEchoLevel(5)
        model_part.CloneTimeStep(1)
        solver.Solve()

        return node_1, node_2, curve

    def testOneLinearSpanWithoutParameterDistortion(self):
        def create_geometry(model_part):
            node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            node_2 = model_part.CreateNewNode(2, 2.0, 0.0, 0.0)

            nodes = KM.NodesVector()
            nodes.append(node_1)
            nodes.append(node_2)

            knots_u = KM.Vector(2)
            knots_u[0] = 0.0
            knots_u[1] = 2.0

            curve = KM.NurbsCurveGeometry3D(nodes, 1, knots_u)

            return node_1, node_2, curve

        node_1, node_2, _ = TrussElementTests.solve(create_geometry)

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

            nodes = KM.NodesVector()
            nodes.append(node_1)
            nodes.append(node_2)

            knots_u = KM.Vector(2)
            knots_u[0] = 0.0
            knots_u[1] = 5.0

            curve = KM.NurbsCurveGeometry3D(nodes, 1, knots_u)

            return node_1, node_2, curve

        node_1, node_2, _ = TrussElementTests.solve(create_geometry)

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

            nodes = KM.NodesVector()
            nodes.append(node_1)
            nodes.append(node_2)
            nodes.append(node_3)

            knots_u = KM.Vector(3)
            knots_u[0] = 0.0
            knots_u[1] = 1.0
            knots_u[2] = 2.0

            curve = KM.NurbsCurveGeometry3D(nodes, 1, knots_u)

            return node_1, node_3, curve

        node_1, node_2, _ = TrussElementTests.solve(create_geometry)

        self.assertAlmostEqual(node_1.X, 0.0)
        self.assertAlmostEqual(node_1.Y, 0.0)
        self.assertAlmostEqual(node_1.Z, 0.0)

        self.assertAlmostEqual(node_2.X, 2.6268671884353934)
        self.assertAlmostEqual(node_2.Y, 0.0)
        self.assertAlmostEqual(node_2.Z, 0.0)

    def testTwoLinearSpansWithParameterDistortion(self):
        def create_geometry(model_part):
            node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
            node_3 = model_part.CreateNewNode(3, 2.0, 0.0, 0.0)

            nodes = KM.NodesVector()
            nodes.append(node_1)
            nodes.append(node_2)
            nodes.append(node_3)

            knots_u = KM.Vector(3)
            knots_u[0] = 0.0
            knots_u[1] = 3.0
            knots_u[2] = 4.0

            curve = KM.NurbsCurveGeometry3D(nodes, 1, knots_u)

            return node_1, node_3, curve

        node_1, node_2, _ = TrussElementTests.solve(create_geometry)

        self.assertAlmostEqual(node_1.X, 0.0)
        self.assertAlmostEqual(node_1.Y, 0.0)
        self.assertAlmostEqual(node_1.Z, 0.0)

        self.assertAlmostEqual(node_2.X, 2.6268671884353934)
        self.assertAlmostEqual(node_2.Y, 0.0)
        self.assertAlmostEqual(node_2.Z, 0.0)

    def testOneQuadraticSpan(self):
        def create_geometry(model_part):
            node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
            node_3 = model_part.CreateNewNode(3, 2.0, 0.0, 0.0)

            nodes = KM.NodesVector()
            nodes.append(node_1)
            nodes.append(node_2)
            nodes.append(node_3)

            knots_u = KM.Vector(4)
            knots_u[0] = 0.0
            knots_u[1] = 0.0
            knots_u[2] = 2.0
            knots_u[3] = 2.0

            curve = KM.NurbsCurveGeometry3D(nodes, 2, knots_u)

            return node_1, node_3, curve

        node_1, node_2, _ = TrussElementTests.solve(create_geometry)

        self.assertAlmostEqual(node_1.X, 0.0)
        self.assertAlmostEqual(node_1.Y, 0.0)
        self.assertAlmostEqual(node_1.Z, 0.0)

        self.assertAlmostEqual(node_2.X, 2.6268671884353934)
        self.assertAlmostEqual(node_2.Y, 0.0)
        self.assertAlmostEqual(node_2.Z, 0.0)
