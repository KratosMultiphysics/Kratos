import KratosMultiphysics as KM
from KratosMultiphysics.IgaApplication import *

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest

def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

class TestNurbsVolumeElement(KratosUnittest.TestCase):
    def create_geometry(model_part, order_u, order_v_w, length):
        nodes = KM.NodesVector()
        # Create control points
        node_id = 1
        width = 1.0
        delta_u = length / order_u
        delta_v = width / order_v_w
        delta_w = width / order_v_w
        z_list = [x for x in drange(0.0, width + 0.5*delta_w, delta_w)]
        y_list = [x for x in drange(0.0, width + 0.5*delta_v, delta_v)]
        x_list = [x for x in drange(0.0, length + 0.5*delta_u, delta_u)]
        for k in z_list: #z-direction
            for j in y_list: #y-direction
                for i in x_list: #x-direction
                    node = model_part.CreateNewNode(node_id, i, j, k)
                    nodes.append(node)
                    node_id = node_id + 1

        knots_u = KM.Vector(2*order_u)
        for i in range(2*order_u):
            if( i < order_u):
                knots_u[i] = 0.0
            else:
                knots_u[i] = 1.0

        knots_v = KM.Vector(2*order_v_w)
        knots_w = KM.Vector(2*order_v_w)
        for i in range(2*order_v_w):
            if( i < order_v_w):
                knots_v[i] = 0.0
                knots_w[i] = 0.0
            else:
                knots_v[i] = 1.0
                knots_w[i] = 1.0

        volume = KM.NurbsVolumeGeometry(
            nodes, order_u, order_v_w, order_v_w, knots_u, knots_v, knots_w)

        return volume

    def solve_cantilever(order_u, order_v_w, length):
        model = KM.Model()
        model_part = model.CreateModelPart('Model')

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)
        model_part.AddNodalSolutionStepVariable(POINT_LOAD)
        model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE,3)

        # Create property for volume elements
        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        volume_properties = model_part.GetProperties()[0]

        volume_properties.SetValue(KM.YOUNG_MODULUS, 210000000)
        volume_properties.SetValue(KM.POISSON_RATIO, 0.3)
        volume_properties.SetValue(KM.DENSITY, 78.5)
        volume_properties.SetValue(KM.CONSTITUTIVE_LAW, cl)

        # Create a nurbs volume
        volume = TestNurbsVolumeElement.create_geometry(model_part, order_u, order_v_w, length)

        # Create quadrature_point_geometries
        quadrature_point_geometries = KM.GeometriesVector()

        volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2)

        id = 1
        for i in range(0, len(quadrature_point_geometries)):
            model_part.CreateNewElement('UpdatedLagrangianElement3D4N', id, quadrature_point_geometries[i], volume_properties)
            id += 1

        # Add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

        for node in volume:
            if( node.X0 == 0.0 ):
                volume[node.Id-1].Fix(KM.DISPLACEMENT_X)
                volume[node.Id-1].Fix(KM.DISPLACEMENT_Y)
                volume[node.Id-1].Fix(KM.DISPLACEMENT_Z)

        # Apply neumann conditions
        prop = model_part.GetProperties()[2]
        force = -1000
        node_id_force = []
        for node in volume:
            if( node.X0 == length):
                node_id_force.append(node.Id)
        nodal_force = force/len(node_id_force)

        id_counter = 1
        for i in node_id_force:
            model_part.CreateNewCondition('PointLoadCondition3D1N', id_counter, [i], prop)
            volume[i-1].SetSolutionStepValue(POINT_LOAD_Y, nodal_force)
            id_counter = id_counter + 1

        # Setup solver
        model_part.SetBufferSize(2)
        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        linear_solver = KM.SkylineLUFactorizationSolver()

        builder_and_solver = KM.ResidualBasedBlockBuilderAndSolver(linear_solver)
        builder_and_solver.SetEchoLevel(0)
        relative_tolerance = 1e-4
        absolute_tolerance = 1e-9
        conv_criteria = KM.ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(0)

        maximum_iterations = 5
        compute_reactions = False
        reform_dofs_at_each_iteration = False
        move_mesh_flag = True

        solver = KM.ResidualBasedNewtonRaphsonStrategy(
            model_part,
            time_scheme,
            conv_criteria,
            builder_and_solver,
            maximum_iterations,
            compute_reactions,
            reform_dofs_at_each_iteration,
            move_mesh_flag
        )

        solver.SetEchoLevel(0)
        model_part.CloneTimeStep(1)
        solver.Solve()

        return volume

    def test_cantilever_beam(self):
        polynomial_degree_u = 4
        polynomial_degree_v_w = 2
        length = 5.0
        volume = TestNurbsVolumeElement.solve_cantilever(polynomial_degree_u, polynomial_degree_v_w, length)

        # Check control points
        self.assertAlmostEqual(volume[4].GetSolutionStepValue(KM.DISPLACEMENT_X),-0.00034683918083963053)
        self.assertAlmostEqual(volume[4].GetSolutionStepValue(KM.DISPLACEMENT_Y),-0.0023516214172830652)
        self.assertAlmostEqual(volume[4].GetSolutionStepValue(KM.DISPLACEMENT_Z),1.8440276704657843e-06)

        self.assertAlmostEqual(volume[9].GetSolutionStepValue(KM.DISPLACEMENT_X),-6.783573692241666e-07)
        self.assertAlmostEqual(volume[9].GetSolutionStepValue(KM.DISPLACEMENT_Y),-0.002355935614369094)
        self.assertAlmostEqual(volume[9].GetSolutionStepValue(KM.DISPLACEMENT_Z),-6.987748079287846e-10)

        self.assertAlmostEqual(volume[14].GetSolutionStepValue(KM.DISPLACEMENT_X),0.0003454751294881693)
        self.assertAlmostEqual(volume[14].GetSolutionStepValue(KM.DISPLACEMENT_Y),-0.002351852289133881)
        self.assertAlmostEqual(volume[14].GetSolutionStepValue(KM.DISPLACEMENT_Z),-1.841706857277202e-06)

        self.assertAlmostEqual(volume[19].GetSolutionStepValue(KM.DISPLACEMENT_X),-0.000352774752785626)
        self.assertAlmostEqual(volume[19].GetSolutionStepValue(KM.DISPLACEMENT_Y),-0.0023518504507103045)
        self.assertAlmostEqual(volume[19].GetSolutionStepValue(KM.DISPLACEMENT_Z),-2.9875801886591146e-15)

        self.assertAlmostEqual(volume[24].GetSolutionStepValue(KM.DISPLACEMENT_X),-6.713337797836725e-07)
        self.assertAlmostEqual(volume[24].GetSolutionStepValue(KM.DISPLACEMENT_Y),-0.0023534019400475166)
        self.assertAlmostEqual(volume[24].GetSolutionStepValue(KM.DISPLACEMENT_Z),-3.0701249540341e-15)

        self.assertAlmostEqual(volume[29].GetSolutionStepValue(KM.DISPLACEMENT_X),0.0003514042028242756)
        self.assertAlmostEqual(volume[29].GetSolutionStepValue(KM.DISPLACEMENT_Y),-0.00235209100722776)
        self.assertAlmostEqual(volume[29].GetSolutionStepValue(KM.DISPLACEMENT_Z),-2.409258133072183e-15)

        self.assertAlmostEqual(volume[34].GetSolutionStepValue(KM.DISPLACEMENT_X),-0.0003468391808390897)
        self.assertAlmostEqual(volume[34].GetSolutionStepValue(KM.DISPLACEMENT_Y),-0.0023516214172832313)
        self.assertAlmostEqual(volume[34].GetSolutionStepValue(KM.DISPLACEMENT_Z),-1.844027676601746e-06)

        self.assertAlmostEqual(volume[39].GetSolutionStepValue(KM.DISPLACEMENT_X),-6.783573681337827e-07)
        self.assertAlmostEqual(volume[39].GetSolutionStepValue(KM.DISPLACEMENT_Y),-0.002355935614369571)
        self.assertAlmostEqual(volume[39].GetSolutionStepValue(KM.DISPLACEMENT_Z),6.987692955943273e-10)

        self.assertAlmostEqual(volume[44].GetSolutionStepValue(KM.DISPLACEMENT_X),0.0003454751294886004)
        self.assertAlmostEqual(volume[44].GetSolutionStepValue(KM.DISPLACEMENT_Y),-0.0023518522891342523)
        self.assertAlmostEqual(volume[44].GetSolutionStepValue(KM.DISPLACEMENT_Z),1.8417068522526938e-06)

        # Check global coordinate
        global_coord = KM.Vector(3)
        param = KM.Vector(3)
        param[0] = 1.0
        param[1] = 0.5
        param[2] = 0.5
        global_coord = volume.GlobalCoordinates(global_coord,param)
        self.assertAlmostEqual(global_coord[0],5.0,5)
        self.assertAlmostEqual(global_coord[1],0.4976467387158514)
        self.assertAlmostEqual(global_coord[2],0.5,5)

if __name__ == '__main__':
    KratosUnittest.main()