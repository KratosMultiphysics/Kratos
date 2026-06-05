#   KRATOS  _____________
#          /  _/ ____/   |
#          / // / __/ /| |
#        _/ // /_/ / ___ |
#       /___/\____/_/  |_| Application

import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest
import os


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class Shell6pElementTests(KratosUnittest.TestCase):

    def solve_cantilever(create_geometry, iterations):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

        # create property for shell elements

        shell_properties = model_part.GetProperties()[1]
        shell_properties.SetValue(KM.THICKNESS, 0.1)
        shell_properties.SetValue(KM.YOUNG_MODULUS, 200000000)
        shell_properties.SetValue(KM.POISSON_RATIO, 0)
        shell_properties.SetValue(KM.CONSTITUTIVE_LAW, SMA.LinearElastic3DLaw())

        # create a nurbs surface
        surface = create_geometry(model_part)

        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()

        surface.CreateQuadraturePointGeometries(quadrature_point_geometries, 3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('Shell6pElement', element_id, quadrature_point_geometries[i], shell_properties)
            element_id += 1

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)

        # apply dirichlet conditions
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

        maximum_iterations = iterations 
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
    
    def solve_cantilever_weak_support(create_geometry, iterations):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

        # create property for shell elements and penalty conditions
        shell_properties = model_part.GetProperties()[1]
        shell_properties.SetValue(KM.THICKNESS, 0.1)
        shell_properties.SetValue(KM.YOUNG_MODULUS, 200000000)
        shell_properties.SetValue(KM.POISSON_RATIO, 0)
        shell_properties.SetValue(KM.CONSTITUTIVE_LAW, SMA.LinearElastic3DLaw())

        penalty_properties = model_part.GetProperties()[2]
        penalty_properties.SetValue(IGA.PENALTY_FACTOR, 1000000000.0)
        penalty_properties.SetValue(IGA.PENALTY_ROTATION_FACTOR, 2000000000.0)
        penalty_properties.SetValue(IGA.INTEGRATE_CONSERVATIVE, True)

        # create a brep surface
        brep_surface, brep_curve_on_surface_left = create_geometry(model_part)

        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()
        quadrature_point_curve_on_surface_geometries = KM.GeometriesVector()

        brep_surface.CreateQuadraturePointGeometries(quadrature_point_geometries, 3)
        surface = brep_surface.GetGeometryPart(KM.Geometry.BACKGROUND_GEOMETRY_INDEX)
    
        brep_curve_on_surface_left.CreateQuadraturePointGeometries(quadrature_point_curve_on_surface_geometries, 3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('Shell6pElement', element_id, quadrature_point_geometries[i], shell_properties)
            element_id += 1

        condition_id = element_id + 1
        for i in range(0, len(quadrature_point_curve_on_surface_geometries)):
            model_part.CreateNewCondition('SupportPenalty6pCondition', condition_id, quadrature_point_curve_on_surface_geometries[i], penalty_properties)
            condition_id += 1

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)

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

        maximum_iterations = iterations 
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

        return brep_surface
    
    def solve_cantilever_modal_analysis(create_geometry):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        model_part.AddNodalSolutionStepVariable(SMA.POINT_LOAD)

        # create property for shell elements

        shell_properties = model_part.GetProperties()[1]
        shell_properties.SetValue(KM.THICKNESS, 0.1)
        shell_properties.SetValue(KM.YOUNG_MODULUS, 200000000)
        shell_properties.SetValue(KM.POISSON_RATIO, 0)
        shell_properties.SetValue(KM.DENSITY, 7850)
        shell_properties.SetValue(KM.CONSTITUTIVE_LAW, SMA.LinearElastic3DLaw())

        # create a nurbs surface
        surface = create_geometry(model_part)

        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()

        surface.CreateQuadraturePointGeometries(quadrature_point_geometries, 3)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('Shell6pElement', element_id, quadrature_point_geometries[i], shell_properties)
            element_id += 1

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)

        # apply dirichlet conditions
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

        # setup solver
        model_part.SetBufferSize(1)
        eigensolver_settings = KM.Parameters("""
        {
            "max_iteration"         : 1000,
            "tolerance"             : 1e-6,
            "number_of_eigenvalues" : 10,
            "echo_level"            : 0,
            "normalize_eigenvectors": true
        }
        """)
        
        eigen_solver = LinearSolversApplication.EigensystemSolver(eigensolver_settings)
        builder_and_solver = KM.ResidualBasedBlockBuilderAndSolver(eigen_solver)
        eigen_scheme = SMA.EigensolverDynamicScheme()
        mass_matrix_diagonal_value = 0.0
        stiffness_matrix_diagonal_value = 1.0

        eigen_solver = SMA.EigensolverStrategy(model_part, eigen_scheme, builder_and_solver,
            mass_matrix_diagonal_value,
            stiffness_matrix_diagonal_value)

        eigen_solver.SetEchoLevel(0)
        eigen_solver.Solve()

        eigenvalues = model_part.ProcessInfo[SMA.EIGENVALUE_VECTOR]

        return surface, eigenvalues

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

        # linear analysis
        surface = Shell6pElementTests.solve_cantilever(create_geometry, 1)

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

    # test same geometry but with different knot vector parameters.
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

        # linear analysis
        surface = Shell6pElementTests.solve_cantilever(create_geometry, 1)

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

    # test same geometry with nonlinear analysis
    def testCantileverOneQuadraticSpanNonlinearAnalysis(self):
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

        # nonlinear analysis
        surface = Shell6pElementTests.solve_cantilever(create_geometry, 100)

        for node in surface:
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_X), 0)
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_Y), 0)

        self.assertAlmostEqual(surface[0].Z, 0.0)
        self.assertAlmostEqual(surface[3].Z, 0.0)
        self.assertAlmostEqual(surface[6].Z, 0.0)

        self.assertAlmostEqual(surface[1].Z, -0.000555264278525)
        self.assertAlmostEqual(surface[4].Z, -0.000447785804707)
        self.assertAlmostEqual(surface[7].Z, -0.000555264278525)

        self.assertAlmostEqual(surface[2].Z, -0.224868627203197)
        self.assertAlmostEqual(surface[5].Z, -0.225040621609457)
        self.assertAlmostEqual(surface[8].Z, -0.224868627203197)

    # test same geometry with modal analysis
    def testCantileverOneQuadraticSpanModalAnalysis(self):
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

        # modal analysis
        surface, eigenvalues = Shell6pElementTests.solve_cantilever_modal_analysis(create_geometry)


    # test nonlinear analysis with the same geometry and weak penalty support
    def testCantileverOneQuadraticSpanWithWeakPenaltySupport(self):
        
        def CreateLinearParametricCurveOnUnitSquare(p0, p1):
            knot_vector_curve = KM.Vector(2)
            knot_vector_curve[0] = 0.0
            knot_vector_curve[1] = 1.0

            degree_curve = 1
            points_curve = MakePointsList2D(p0, p1)

            return KM.NurbsCurveGeometry2DPoint(
                points_curve,
                degree_curve,
                knot_vector_curve
            )

        def MakeBoundaryBrepCurve(surface, p0, p1, same_curve_direction=True):
            curve_2d = CreateLinearParametricCurveOnUnitSquare(p0, p1)
            brep_curve = KM.BrepCurveOnSurface(surface, curve_2d, same_curve_direction)
            # brep_curve._curve_2d = curve_2d
            return brep_curve
        
        def MakePointsList2D(p0, p1):
            return [
                KM.Point(float(p0[0]), float(p0[1]), 0.0),
                KM.Point(float(p1[0]), float(p1[1]), 0.0)
            ]
    
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
            surface.SetId(1001)

            u_min = float(knots_u[0])
            u_max = float(knots_u[len(knots_u) - 1])
            v_min = float(knots_v[0])
            v_max = float(knots_v[len(knots_v) - 1])

            curve_bottom = MakeBoundaryBrepCurve(surface, (u_min, v_min), (u_max, v_min), True)
            curve_right  = MakeBoundaryBrepCurve(surface, (u_max, v_min), (u_max, v_max), True)
            curve_top    = MakeBoundaryBrepCurve(surface, (u_max, v_max), (u_min, v_max), True)
            curve_left   = MakeBoundaryBrepCurve(surface, (u_min, v_max), (u_min, v_min), True)

            outer_loops = [[curve_bottom, curve_right, curve_top, curve_left]]
            inner_loops = []

            brep_surface = KM.BrepSurface(surface, outer_loops, inner_loops)
            brep_surface.SetId(1002) 

            brep_curve_on_surface_left = MakeBoundaryBrepCurve(surface, (u_min, v_max), (u_min, v_min), True)
            brep_curve_on_surface_left.SetId(1003) 
            
            return brep_surface, brep_curve_on_surface_left

        # nonlinear analysis
        brep_surface = Shell6pElementTests.solve_cantilever_weak_support(create_geometry, 100)
        surface = brep_surface.GetGeometryPart(KM.Geometry.BACKGROUND_GEOMETRY_INDEX)

        for node in surface:
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_X), 0)
            self.assertAlmostEqual(node.GetValue(KM.DISPLACEMENT_Y), 0)

        self.assertAlmostEqual(surface[0].Z, -1.22636367511e-07)
        self.assertAlmostEqual(surface[3].Z, -3.54727264976e-07)
        self.assertAlmostEqual(surface[6].Z, -1.22636367510e-07)

        self.assertAlmostEqual(surface[1].Z, -0.000556853191498)
        self.assertAlmostEqual(surface[4].Z, -0.000448953526144)
        self.assertAlmostEqual(surface[7].Z, -0.000556853191498)

        self.assertAlmostEqual(surface[2].Z, -0.224871260111161)
        self.assertAlmostEqual(surface[5].Z, -0.225043464084019)
        self.assertAlmostEqual(surface[8].Z, -0.224871260111161)