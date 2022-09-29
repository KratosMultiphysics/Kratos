import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.gid_output_process import GiDOutputProcess

if kratos_utilities.CheckIfApplicationsAvailable("LinearSolversApplication"):
    import KratosMultiphysics.LinearSolversApplication as KratosLinearSolvers

class TestRBFShapeFunctionsUtility(KratosUnittest.TestCase):

    def setUp(self):
        # Create the test model part
        self.model = KratosMultiphysics.Model()
        self.main_model_part = self.model.CreateModelPart("MainModelPart")
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        # Generate the problem domain
        problem_domain = KratosMultiphysics.Quadrilateral2D4(
            KratosMultiphysics.Node(1, 0.0, 0.0, 0.0),
            KratosMultiphysics.Node(2, 0.0, 1.0, 0.0),
            KratosMultiphysics.Node(3, 1.0, 1.0, 0.0),
            KratosMultiphysics.Node(4, 1.0, 0.0, 0.0))

        mesh_parameters = KratosMultiphysics.Parameters("{}")
        mesh_parameters.AddEmptyValue("element_name").SetString("Element2D3N")
        mesh_parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
        mesh_parameters.AddEmptyValue("number_of_divisions").SetInt(5)

        KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, self.main_model_part, mesh_parameters).Execute()

    def test_rbf_shape_functions_utility_1x1_square(self):
        # Create the sample points array
        n_nodes = self.main_model_part.NumberOfNodes()
        pts_coord = KratosMultiphysics.Matrix(n_nodes, 3, 0.0)
        i = 0
        for node in self.main_model_part.Nodes:
            pts_coord[i,0] = node.X
            pts_coord[i,1] = node.Y
            i += 1

        # Calculate the RBF shape functions in the square midpoint
        midpoint = KratosMultiphysics.Vector(3)
        midpoint[0] = 0.5
        midpoint[1] = 0.5
        midpoint[2] = 0.0
        N_container = KratosMultiphysics.Vector()
        h = 1 #Shape parameter for the kernel function (radial basis function, inverse multiquadratic)
        input_shape_parameter = False
        if input_shape_parameter:
            KratosMultiphysics.RBFShapeFunctionsUtility.CalculateShapeFunctions(
                pts_coord,
                midpoint,
                h,
                N_container)
        else:
            KratosMultiphysics.RBFShapeFunctionsUtility.CalculateShapeFunctions(
                pts_coord,
                midpoint,
                N_container)

        '''
        For debugging:
        # Save the obtained results in TEMPERATURE variable and visualize
        output_results = False
        if output_results:
            i = 0
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, N_container[i])
                i += 1
            self.PostProcess("test_rbf_shape_functions_utility_1x1_square")
        '''

        # Check results
        N_tot = sum(N_container)
        if input_shape_parameter:
            self.assertAlmostEqual(N_tot, 1.0, 8)
            self.assertAlmostEqual(N_container[0], -0.00044077745800828917, 8)
            self.assertAlmostEqual(N_container[9], -0.08250832642936118, 8)
            self.assertAlmostEqual(N_container[15], 0.3675458064498694, 8)
        else:
            self.assertAlmostEqual(N_tot, 1.0, 8)
            self.assertAlmostEqual(N_container[0], -0.000584922796671384, 8)
            self.assertAlmostEqual(N_container[9], -0.04990247800841474, 8)
            self.assertAlmostEqual(N_container[15], 0.332227735376926, 8)

    @KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
    def test_rbf_shape_functions_utility_ill_conditioned(self):
        # Create the sample points array
        # Note that we deliberately duplicate one point to make the case ill-conditioned
        cloud_points_array = [
            [1.0,0.0,0.0],
            [0.0,1.0,0.0],
            [-1.0,0.0,0.0],
            [-0.99,0.0,0.0]]
        cloud_coords = KratosMultiphysics.Matrix(len(cloud_points_array), 3)
        for i in range(len(cloud_points_array)):
            for j in range(3):
                cloud_coords[i, j] = cloud_points_array[i][j]

        # Calculate the RBF shape functions in the coordinates origin
        # Note that we use the QR with pivoting from Eigen in order to solve the ill-conditioned system
        N_container = KratosMultiphysics.Vector()
        point_coords = KratosMultiphysics.Vector(3, 0.0)
        eigen_col_piv_qr = KratosLinearSolvers.EigenDenseColumnPivotingHouseholderQRDecomposition()
        KratosMultiphysics.RBFShapeFunctionsUtility.CalculateShapeFunctions(
            cloud_coords,
            point_coords,
            N_container,
            eigen_col_piv_qr)

        N_expected = KratosMultiphysics.Vector(4)
        N_expected[0] = 0.3477832764687009
        N_expected[1] = 0.2813057493661467
        N_expected[2] = -11.243513620898828
        N_expected[3] = 11.61442459506398
        self.assertVectorAlmostEqual(N_container, N_expected, 8)

    def PostProcess(self, output_name):
        gid_output = GiDOutputProcess(
            self.main_model_part,
            output_name,
            KratosMultiphysics.Parameters("""{
                "result_file_configuration" : {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostBinary",
                        "WriteDeformedMeshFlag": "WriteUndeformed",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "nodal_results" : ["TEMPERATURE"]
                }
            }"""))
        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()

if __name__ == '__main__':
    KratosUnittest.main()
