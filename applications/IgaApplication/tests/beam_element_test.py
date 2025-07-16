import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest
import numpy as np

def compute_t0_and_n0(curve):
    t0 = np.array(curve.GlobalSpaceDerivatives(
        [KM.Array3([0, 0, 0])], KM.Array3([0.0, 0.0, 0.0]), 1)[1])
    t0 /= np.linalg.norm(t0)
    n0 = np.cross(np.array([0, -1, 0], dtype=float), t0)
    return n0, t0

def w_analytical(x, F, E, I, L):
    return F / (6 * E * I) * (3 * L * x**2 - x**3)

class Beam4pElementTest(KratosUnittest.TestCase):

    def solve(create_geometry):
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
        beam_properties.SetValue(IGA.I_Y, w * h**3 / 12)
        beam_properties.SetValue(IGA.I_Z, h * w**3 / 12)
        beam_properties.SetValue(KM.YOUNG_MODULUS   , 70000)
        beam_properties.SetValue(KM.POISSON_RATIO   , 0)
        beam_properties.SetValue(KM.DENSITY         , 7856)
        beam_properties.SetValue(IGA.CENTER_LINE_ROTATION, np.array([[0, 0], [0.25 , 0], [0.5 , 0],[0.75 , 0], [1, 0]]))

        # create a node based geometry
        node_1, node_2, last_node, curve = create_geometry(model_part)
        n0, t0 = compute_t0_and_n0(curve)
        beam_properties.SetValue(IGA.T_0, t0)
        beam_properties.SetValue(IGA.N_0, n0)
        
        # create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()
        curve.CreateQuadraturePointGeometries(quadrature_point_geometries, 10)

        element_id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('IsogeometricBeamElement', element_id, quadrature_point_geometries[i], beam_properties)
            element_id += 1

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, model_part)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, model_part)# dummy DOF
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, model_part)# dummy DOF

        # apply dirichlet conditions
        node_1.Fix(KM.DISPLACEMENT_X)
        node_1.Fix(KM.DISPLACEMENT_Y)
        node_1.Fix(KM.DISPLACEMENT_Z)

        node_2.Fix(KM.DISPLACEMENT_X)
        node_2.Fix(KM.DISPLACEMENT_Y)
        node_2.Fix(KM.DISPLACEMENT_Z)

        # apply neumann conditions
        prop = model_part.GetProperties()[2]
        model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [last_node.Id], prop)
        last_node.SetSolutionStepValue(SMA.POINT_LOAD, KM.Array3([0,  0, -1]))

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

        solver.SetEchoLevel(5)
        model_part.CloneTimeStep(1)
        solver.Solve()

        return node_1, node_2, last_node, curve

    def testClampedWithoutParameterDistortion(self):
        def create_geometry(model_part):
            node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
            node_3 = model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
            node_4 = model_part.CreateNewNode(4, 3.0, 0.0, 0.0)
            node_5 = model_part.CreateNewNode(5, 4.0, 0.0, 0.0)

            nodes = KM.NodesVector()
            nodes.append(node_1)
            nodes.append(node_2)
            nodes.append(node_3)
            nodes.append(node_4)
            nodes.append(node_5)

            knots = KM.Vector(8)
            knots[0] = 0
            knots[1] = 0
            knots[2] = 0
            knots[3] = 0
            knots[4] = 1
            knots[5] = 1
            knots[6] = 1
            knots[7] = 1


            curve = KM.NurbsCurveGeometry3D(nodes, 4, knots)
            return node_1, node_2, node_5, curve

        node_1, node_2, last_node , _ = Beam4pElementTest.solve(create_geometry)

        self.assertAlmostEqual(node_1.X, 0.0               )
        self.assertAlmostEqual(node_1.Y, 0.0               )
        self.assertAlmostEqual(node_1.Z, 0.0               )

        self.assertAlmostEqual(node_2.X, 1.0               )
        self.assertAlmostEqual(node_2.Y, 0.0               )
        self.assertAlmostEqual(node_2.Z, 0.0               )

        self.assertAlmostEqual(np.array(last_node.GetSolutionStepValue(KM.DISPLACEMENT))[0], 0.0               )
        self.assertAlmostEqual(np.array(last_node.GetSolutionStepValue(KM.DISPLACEMENT))[1], 0.0               )
        self.assertAlmostEqual(np.array(last_node.GetSolutionStepValue(KM.DISPLACEMENT))[2], w_analytical(last_node.X,-1, 70000,5**3*2.5/12,4))
