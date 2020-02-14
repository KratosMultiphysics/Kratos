#
#   KRATOS  _____________
#          /  _/ ____/   |
#          / // / __/ /| |
#        _/ // /_/ / ___ |
#       /___/\____/_/  |_| Application
#


import KratosMultiphysics as KM
from KratosMultiphysics.IgaApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class Shell3pElementTests(KratosUnittest.TestCase):

    def solve_cantilever(create_geometry):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(POINT_LOAD)

        # create property for shell elements

        shell_properties = model_part.GetProperties()[1]
        shell_properties.SetValue(KM.THICKNESS       , 0.1      )
        shell_properties.SetValue(KM.YOUNG_MODULUS   , 210000000)
        shell_properties.SetValue(KM.POISSON_RATIO   , 0        )
        shell_properties.SetValue(KM.DENSITY         , 78.5     )
        shell_properties.SetValue(KM.CONSTITUTIVE_LAW,
            LinearElasticPlaneStress2DLaw())

        # create a nurbs surface
        surface = create_geometry(model_part)

        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()

        surface.CreateQuadraturePointGeometries(quadrature_point_geometries, 3)

        KM.CreateElements(quadrature_point_geometries, model_part, 'Shell3pElement', 1, 0)

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

        # apply dirichlet conditions
        surface[0].Fix(KM.DISPLACEMENT_X)
        surface[0].Fix(KM.DISPLACEMENT_Y)
        surface[0].Fix(KM.DISPLACEMENT_Z)

        surface[1].Fix(KM.DISPLACEMENT_X)
        surface[1].Fix(KM.DISPLACEMENT_Y)
        surface[1].Fix(KM.DISPLACEMENT_Z)

        surface[2].Fix(KM.DISPLACEMENT_X)
        surface[2].Fix(KM.DISPLACEMENT_Y)
        surface[2].Fix(KM.DISPLACEMENT_Z)

        # clamp bending
        surface[3].Fix(KM.DISPLACEMENT_Z)
        surface[4].Fix(KM.DISPLACEMENT_Z)
        surface[5].Fix(KM.DISPLACEMENT_Z)

        # apply neumann conditions

        prop = model_part.GetProperties()[2]

        model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [7], prop)
        surface[6].SetSolutionStepValue(POINT_LOAD_Z, -50)

        model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [9], prop)
        surface[8].SetSolutionStepValue(POINT_LOAD_Z, -50)

        # setup solver
        model_part.SetBufferSize(1)
        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        linear_solver = linear_solver_factory.ConstructSolver(KM.Parameters(
            r'{"solver_type": "skyline_lu_factorization"}'))

        relative_tolerance = 1e-8
        absolute_tolerance = 1e-7

        conv_criteria = KM.ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(0)

        maximum_iterations = 100
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

        surface = Shell3pElementTests.solve_cantilever(create_geometry)
        
        for node in surface:
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_X), 0)
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_Y), 0)

        self.assertAlmostEqual(surface[0].Z, 0.0)
        self.assertAlmostEqual(surface[1].Z, 0.0)
        self.assertAlmostEqual(surface[2].Z, 0.0)

        self.assertAlmostEqual(surface[3].Z, 0.0)
        self.assertAlmostEqual(surface[4].Z, 0.0)
        self.assertAlmostEqual(surface[5].Z, 0.0)

        self.assertAlmostEqual(surface[6].Z, -0.141013552223511)
        self.assertAlmostEqual(surface[7].Z, -0.140899868451905)
        self.assertAlmostEqual(surface[8].Z, -0.141013552223511)

    # test same geometry but with different knot  vector parameters.
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

        surface = ShellKLDiscreteElementTests.solve_cantilever(create_geometry)

        for node in surface:
            self.assertAlmostEqual(node.GetValue(DISPLACEMENT_X), 0)
            self.assertAlmostEqual(node.GetValue(DISPLACEMENT_Y), 0)

        self.assertAlmostEqual(surface[0].Z, 0.0)
        self.assertAlmostEqual(surface[1].Z, 0.0)
        self.assertAlmostEqual(surface[2].Z, 0.0)

        self.assertAlmostEqual(surface[3].Z, 0.0)
        self.assertAlmostEqual(surface[4].Z, 0.0)
        self.assertAlmostEqual(surface[5].Z, 0.0)

        self.assertAlmostEqual(surface[6].Z, -0.141013552223511)
        self.assertAlmostEqual(surface[7].Z, -0.140899868451905)
        self.assertAlmostEqual(surface[8].Z, -0.141013552223511)
