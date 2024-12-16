import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.gid_output_process import GiDOutputProcess

def PostProcess(model_part, output_name):
    gid_output = GiDOutputProcess(
        model_part,
        output_name,
        KratosMultiphysics.Parameters("""{
            "result_file_configuration" : {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostBinary",
                    "WriteDeformedMeshFlag": "WriteUndeformed",
                    "WriteConditionsFlag": "WriteConditions",
                    "MultiFileFlag": "SingleFile"
                },
                "nodal_results" : ["TEMPERATURE","DISPLACEMENT"]
            }
        }"""))
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

class TestMLSShapeFunctionsUtility(KratosUnittest.TestCase):

    def setUp(self):
        # Create the test model part
        self.model = KratosMultiphysics.Model()
        self.main_model_part = self.model.CreateModelPart("MainModelPart")
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # Generate the problem domain
        problem_domain = KratosMultiphysics.Quadrilateral2D4(
            KratosMultiphysics.Node(1, 0.0, 0.0, 0.0),
            KratosMultiphysics.Node(2, 0.0, 1.0, 0.0),
            KratosMultiphysics.Node(3, 1.0, 1.0, 0.0),
            KratosMultiphysics.Node(4, 1.0, 0.0, 0.0))

        mesh_parameters = KratosMultiphysics.Parameters("{}")
        mesh_parameters.AddEmptyValue("element_name").SetString("Element2D3N")
        mesh_parameters.AddEmptyValue("condition_name").SetString("LineCondition2D2N")
        mesh_parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
        mesh_parameters.AddEmptyValue("number_of_divisions").SetInt(20)

        KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, self.main_model_part, mesh_parameters).Execute()

    def testMLSShapeFunctionsUtility1x1Square1stOrder(self):
        # Create the sample points array
        n_nodes = self.main_model_part.NumberOfNodes()
        pts_coord = KratosMultiphysics.Matrix(n_nodes, 3, 0.0)
        i = 0
        for node in self.main_model_part.Nodes:
            pts_coord[i,0] = node.X
            pts_coord[i,1] = node.Y
            i += 1

        # Calculate the MLS shape functions in the square midpoint
        h = 0.2
        dim = 2
        order = 1
        midpoint = KratosMultiphysics.Vector(3)
        midpoint[0] = 0.5
        midpoint[1] = 0.5
        midpoint[2] = 0.0
        N_container = KratosMultiphysics.Vector()
        DN_DX_container = KratosMultiphysics.Matrix()
        KratosMultiphysics.MLSShapeFunctionsUtility.CalculateShapeFunctionsAndGradients(
            dim,
            order,
            pts_coord,
            midpoint,
            h,
            N_container,
            DN_DX_container)

        # Save the obtained results in TEMPERATURE variable and visualize
        output_results = False
        if output_results:
            i = 0
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, N_container[i])
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, DN_DX_container[i,0])
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, DN_DX_container[i,1])
                i += 1
            PostProcess(self.main_model_part, "testMLSShapeFunctionsUtility1x1Square1stOrder")

        # Check results
        N_tot = sum(N_container)
        self.assertAlmostEqual(N_tot, 1.0, 8)
        self.assertAlmostEqual(N_container[0], 7.416763267826954e-08, 8)
        self.assertAlmostEqual(DN_DX_container[0,0], 4.6903874434833715e-06, 8)
        self.assertAlmostEqual(DN_DX_container[0,1], 4.690387443475345e-06, 8)
        self.assertAlmostEqual(N_container[220], 0.019901941300573312, 8)
        self.assertAlmostEqual(DN_DX_container[220,0], 0.0, 8)
        self.assertAlmostEqual(DN_DX_container[220,1], 0.0, 8)
        self.assertAlmostEqual(N_container[221], 0.01869614363310557, 8)
        self.assertAlmostEqual(DN_DX_container[221,0], 0.0, 8)
        self.assertAlmostEqual(DN_DX_container[221,1], 0.24155121872335386, 8)

    def testMLSShapeFunctionsUtility1x1Square2ndOrder(self):
        # Create the sample points array
        n_nodes = self.main_model_part.NumberOfNodes()
        pts_coord = KratosMultiphysics.Matrix(n_nodes, 3, 0.0)
        i = 0
        for node in self.main_model_part.Nodes:
            pts_coord[i,0] = node.X
            pts_coord[i,1] = node.Y
            i += 1

        # Calculate the MLS shape functions in the square midpoint
        h = 0.2
        dim = 2
        order = 2
        midpoint = KratosMultiphysics.Vector(3)
        midpoint[0] = 0.5
        midpoint[1] = 0.5
        midpoint[2] = 0.0
        N_container = KratosMultiphysics.Vector()
        DN_DX_container = KratosMultiphysics.Matrix()
        KratosMultiphysics.MLSShapeFunctionsUtility.CalculateShapeFunctionsAndGradients(
            dim,
            order,
            pts_coord,
            midpoint,
            h,
            N_container,
            DN_DX_container)

        # Save the obtained results in TEMPERATURE variable and visualize
        output_results = False
        if output_results:
            i = 0
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, N_container[i])
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, DN_DX_container[i,0])
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, DN_DX_container[i,1])
                i += 1
            PostProcess(self.main_model_part, "testMLSShapeFunctionsUtility1x1Square2ndOrder")

        # Check results
        N_tot = sum(N_container)
        self.assertAlmostEqual(N_tot, 1.0, 8)
        self.assertAlmostEqual(N_container[0], -7.951619456545074e-07, 8)
        self.assertAlmostEqual(DN_DX_container[0,0], 2.775864113372096e-06, 8)
        self.assertAlmostEqual(DN_DX_container[0,1], 2.775864113372096e-06, 8)
        self.assertAlmostEqual(N_container[220], 0.04012381031655152, 8)
        self.assertAlmostEqual(DN_DX_container[220,0], 0.0, 8)
        self.assertAlmostEqual(DN_DX_container[220,1], 0.0, 8)
        self.assertAlmostEqual(N_container[221], 0.036502147047586114, 8)
        self.assertAlmostEqual(DN_DX_container[221,0], 0.0, 8)
        self.assertAlmostEqual(DN_DX_container[221,1], 0.6415689416490207, 8)

if __name__ == '__main__':
    KratosUnittest.main()

