#
#   KRATOS  _____________
#          /  _/ ____/   |
#          / // / __/ /| |
#        _/ // /_/ / ___ |
#       /___/\____/_/  |_| Application
#
#   Main authors:   Thomas Oberbichler
#

from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import new_linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class ShellKLDiscreteElementTest(KratosUnittest.TestCase):

    def solve_cantilever(create_geometry):
        model_part = ModelPart('Model')

        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(REACTION)
        model_part.AddNodalSolutionStepVariable(POINT_LOAD)

        # create property for shell elements

        shell_properties = model_part.GetProperties()[1]
        shell_properties.SetValue(THICKNESS       , 0.1      )
        shell_properties.SetValue(YOUNG_MODULUS   , 210000000)
        shell_properties.SetValue(POISSON_RATIO   , 0        )
        shell_properties.SetValue(DENSITY         , 78.5     )
        shell_properties.SetValue(CONSTITUTIVE_LAW,
            LinearElasticPlaneStress2DLaw())

        # create a node based geometry

        surface = create_geometry(model_part)

        # create elements for each integration point

        spans_u = surface.SpansU()
        spans_v = surface.SpansV()

        integration_points = [integration_point
            for integration_point in IntegrationPoints.Points2D(
            surface.DegreeU + 1, surface.DegreeV + 1, spans_u, spans_v)]

        shapes = SurfaceShapeEvaluator(DegreeU=surface.DegreeU,
            DegreeV=surface.DegreeV, Order=2)

        for i, (u, v, weight) in enumerate(integration_points):
            shapes.Compute(surface.KnotsU(), surface.KnotsV(), u, v)

            node_indices = [surface.Node(u, v).Id
                for u, v in shapes.NonzeroPoleIndices]

            element = model_part.CreateNewElement('ShellKLDiscreteElement',
                i + 1, node_indices, shell_properties)

            n_0 = Vector(shapes.NumberOfNonzeroPoles)
            n_1 = Matrix(shapes.NumberOfNonzeroPoles, 2)
            n_2 = Matrix(shapes.NumberOfNonzeroPoles, 3)

            for j, (poleU, poleV) in enumerate(shapes.NonzeroPoleIndices):
                indexU = poleU - shapes.FirstNonzeroPoleU
                indexV = poleV - shapes.FirstNonzeroPoleV
                n_0[j   ] = shapes(0, indexU, indexV)
                n_1[j, 0] = shapes(1, indexU, indexV)
                n_1[j, 1] = shapes(2, indexU, indexV)
                n_2[j, 0] = shapes(3, indexU, indexV)
                n_2[j, 1] = shapes(5, indexU, indexV)
                n_2[j, 2] = shapes(4, indexU, indexV)

            element.SetValue(INTEGRATION_WEIGHT, weight)
            element.SetValue(SHAPE_FUNCTION_VALUES, n_0)
            element.SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, n_1)
            element.SetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES, n_2)

        # add dofs

        VariableUtils().AddDof(DISPLACEMENT_X, REACTION_X, model_part)
        VariableUtils().AddDof(DISPLACEMENT_Y, REACTION_Y, model_part)
        VariableUtils().AddDof(DISPLACEMENT_Z, REACTION_Z, model_part)

        # apply dirichlet conditions

        surface.Node(0, 0).Fix(DISPLACEMENT_X)
        surface.Node(0, 0).Fix(DISPLACEMENT_Y)
        surface.Node(0, 0).Fix(DISPLACEMENT_Z)

        surface.Node(0, 1).Fix(DISPLACEMENT_X)
        surface.Node(0, 1).Fix(DISPLACEMENT_Y)
        surface.Node(0, 1).Fix(DISPLACEMENT_Z)

        surface.Node(0, 2).Fix(DISPLACEMENT_X)
        surface.Node(0, 2).Fix(DISPLACEMENT_Y)
        surface.Node(0, 2).Fix(DISPLACEMENT_Z)

        surface.Node(1, 0).Fix(DISPLACEMENT_Z)

        surface.Node(1, 1).Fix(DISPLACEMENT_Z)

        surface.Node(1, 2).Fix(DISPLACEMENT_Z)

        # apply neumann conditions

        prop = model_part.GetProperties()[2]

        node_n0 = surface.Node(surface.NumberOfPolesU - 1, 0)
        node_nm = surface.Node(surface.NumberOfPolesU - 1, surface.NumberOfPolesV - 1)

        model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [node_n0.Id],
            prop)
        node_n0.SetSolutionStepValue(POINT_LOAD_Z, -50)

        model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [node_nm.Id],
            prop)
        node_nm.SetSolutionStepValue(POINT_LOAD_Z, -50)

        # setup solver

        model_part.SetBufferSize(1)

        time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
            r'{"solver_type": "SkylineLUFactorizationSolver"}'))

        relative_tolerance = 1e-8
        absolute_tolerance = 1e-7

        conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(0)

        maximum_iterations = 100
        compute_reactions = False
        reform_dofs_at_each_iteration = False
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

        return surface

    def testCantileverOneQuadraticSpanWithoutParameterDistortion(self):
        def create_geometry(model_part):
            node_00 = model_part.CreateNewNode(1, 0.0, 0.0, 0)
            node_10 = model_part.CreateNewNode(2, 2.5, 0.0, 0)
            node_20 = model_part.CreateNewNode(3, 5.0, 0.0, 0)

            node_01 = model_part.CreateNewNode(4, 0.0, 0.5, 0)
            node_11 = model_part.CreateNewNode(5, 2.5, 0.5, 0)
            node_21 = model_part.CreateNewNode(6, 5.0, 0.5, 0)

            node_02 = model_part.CreateNewNode(7, 0.0, 1.0, 0)
            node_12 = model_part.CreateNewNode(8, 2.5, 1.0, 0)
            node_22 = model_part.CreateNewNode(9, 5.0, 1.0, 0)

            surface = NodeSurfaceGeometry3D(
                DegreeU=2,
                DegreeV=2,
                NumberOfNodesU=3,
                NumberOfNodesV=3)

            surface.SetKnotU(Index=0, Value=0.0)
            surface.SetKnotU(Index=1, Value=0.0)
            surface.SetKnotU(Index=2, Value=5.0)
            surface.SetKnotU(Index=3, Value=5.0)

            surface.SetKnotV(Index=0, Value=0.0)
            surface.SetKnotV(Index=1, Value=0.0)
            surface.SetKnotV(Index=2, Value=1.0)
            surface.SetKnotV(Index=3, Value=1.0)

            surface.SetNode(IndexU=0, IndexV=0, Value=node_00)
            surface.SetNode(IndexU=1, IndexV=0, Value=node_10)
            surface.SetNode(IndexU=2, IndexV=0, Value=node_20)

            surface.SetNode(IndexU=0, IndexV=1, Value=node_01)
            surface.SetNode(IndexU=1, IndexV=1, Value=node_11)
            surface.SetNode(IndexU=2, IndexV=1, Value=node_21)

            surface.SetNode(IndexU=0, IndexV=2, Value=node_02)
            surface.SetNode(IndexU=1, IndexV=2, Value=node_12)
            surface.SetNode(IndexU=2, IndexV=2, Value=node_22)

            return surface

        surface = ShellKLDiscreteElementTest.solve_cantilever(create_geometry)

        for i in range(surface.NumberOfPolesU):
            for j in range(surface.NumberOfPolesV):
                node = surface.Node(i, j)
                self.assertAlmostEqual(node.GetValue(DISPLACEMENT_X), 0)
                self.assertAlmostEqual(node.GetValue(DISPLACEMENT_Y), 0)

        self.assertAlmostEqual(surface.Node(0, 0).Z, 0.0)
        self.assertAlmostEqual(surface.Node(0, 1).Z, 0.0)
        self.assertAlmostEqual(surface.Node(0, 2).Z, 0.0)

        self.assertAlmostEqual(surface.Node(1, 0).Z, 0.0)
        self.assertAlmostEqual(surface.Node(1, 1).Z, 0.0)
        self.assertAlmostEqual(surface.Node(1, 2).Z, 0.0)

        self.assertAlmostEqual(surface.Node(2, 0).Z, -0.141013552223511)
        self.assertAlmostEqual(surface.Node(2, 1).Z, -0.140899868451905)
        self.assertAlmostEqual(surface.Node(2, 2).Z, -0.141013552223511)

    def testCantileverOneQuadraticSpanWithParameterDistortion(self):
        def create_geometry(model_part):
            node_00 = model_part.CreateNewNode(1, 0.0, 0.0, 0)
            node_10 = model_part.CreateNewNode(2, 2.5, 0.0, 0)
            node_20 = model_part.CreateNewNode(3, 5.0, 0.0, 0)

            node_01 = model_part.CreateNewNode(4, 0.0, 0.5, 0)
            node_11 = model_part.CreateNewNode(5, 2.5, 0.5, 0)
            node_21 = model_part.CreateNewNode(6, 5.0, 0.5, 0)

            node_02 = model_part.CreateNewNode(7, 0.0, 1.0, 0)
            node_12 = model_part.CreateNewNode(8, 2.5, 1.0, 0)
            node_22 = model_part.CreateNewNode(9, 5.0, 1.0, 0)

            surface = NodeSurfaceGeometry3D(
                DegreeU=2,
                DegreeV=2,
                NumberOfNodesU=3,
                NumberOfNodesV=3)

            surface.SetKnotU(Index=0, Value=0.0)
            surface.SetKnotU(Index=1, Value=0.0)
            surface.SetKnotU(Index=2, Value=1.0)
            surface.SetKnotU(Index=3, Value=1.0)

            surface.SetKnotV(Index=0, Value=0.0)
            surface.SetKnotV(Index=1, Value=0.0)
            surface.SetKnotV(Index=2, Value=2.0)
            surface.SetKnotV(Index=3, Value=2.0)

            surface.SetNode(IndexU=0, IndexV=0, Value=node_00)
            surface.SetNode(IndexU=1, IndexV=0, Value=node_10)
            surface.SetNode(IndexU=2, IndexV=0, Value=node_20)

            surface.SetNode(IndexU=0, IndexV=1, Value=node_01)
            surface.SetNode(IndexU=1, IndexV=1, Value=node_11)
            surface.SetNode(IndexU=2, IndexV=1, Value=node_21)

            surface.SetNode(IndexU=0, IndexV=2, Value=node_02)
            surface.SetNode(IndexU=1, IndexV=2, Value=node_12)
            surface.SetNode(IndexU=2, IndexV=2, Value=node_22)

            return surface

        surface = ShellKLDiscreteElementTest.solve_cantilever(create_geometry)

        for i in range(surface.NumberOfPolesU):
            for j in range(surface.NumberOfPolesV):
                node = surface.Node(i, j)
                self.assertAlmostEqual(node.GetValue(DISPLACEMENT_X), 0)
                self.assertAlmostEqual(node.GetValue(DISPLACEMENT_Y), 0)

        self.assertAlmostEqual(surface.Node(0, 0).Z, 0.0)
        self.assertAlmostEqual(surface.Node(0, 1).Z, 0.0)
        self.assertAlmostEqual(surface.Node(0, 2).Z, 0.0)

        self.assertAlmostEqual(surface.Node(1, 0).Z, 0.0)
        self.assertAlmostEqual(surface.Node(1, 1).Z, 0.0)
        self.assertAlmostEqual(surface.Node(1, 2).Z, 0.0)

        self.assertAlmostEqual(surface.Node(2, 0).Z, -0.141013552223511)
        self.assertAlmostEqual(surface.Node(2, 1).Z, -0.140899868451905)
        self.assertAlmostEqual(surface.Node(2, 2).Z, -0.141013552223511)
