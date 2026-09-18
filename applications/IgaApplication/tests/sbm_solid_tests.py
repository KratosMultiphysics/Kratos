import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_creation_utility import TestCreationUtility

def run_modelers(current_model, modelers_list):
    from KratosMultiphysics.modeler_factory import KratosModelerFactory
    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()

class SbmSolidTests(KratosUnittest.TestCase):

    def testSbmSolidCondition(self):
        current_model = KM.Model()
        skin_model_part_outer_initial = current_model.CreateModelPart("skinModelPart_outer_initial")
        skin_model_part_outer_initial.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        skin_model_part_outer_initial.CreateNewProperties(1)
        skin_model_part_outer_initial.CreateNewNode(1, 0.0, 0.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(2, 2.0, 0.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(3, 2.0, 2.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(4, 0.0, 2.0, 0.0)
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_model_part_outer_initial.GetProperties()[1])
        
        skin_model_part = current_model.CreateModelPart("skin_model_part")
        skin_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        
        modeler_settings = KM.Parameters("""
        [
            {
            "modeler_name": "NurbsGeometryModelerSbm",
            "Parameters": {
                    "model_part_name" : "IgaModelPart",
                    "lower_point_xyz": [0.0,0.0,0.0],
                    "upper_point_xyz": [2.0,2.0,0.0],
                    "lower_point_uvw": [0.0,0.0,0.0],
                    "upper_point_uvw": [2.0,2.0,0.0],
                    "polynomial_order" : [3, 3],
                    "number_of_knot_spans" : [5,5],
                    "lambda_outer": 0.5,
                    "number_of_inner_loops": 0,
                    "skin_model_part_outer_initial_name": "skinModelPart_outer_initial",           
                    "skin_model_part_name": "skin_model_part"
                }
            },
            {
                "modeler_name": "IgaModelerSbm",
                "Parameters": {
                    "echo_level": 0,
                    "skin_model_part_name": "skin_model_part",
                    "analysis_model_part_name": "IgaModelPart",
                    "element_condition_list": [
                        {
                            "geometry_type": "SurfaceEdge",
                            "iga_model_part": "SBM_Support_outer",
                            "type": "condition",
                            "name": "SbmSolidCondition",
                            "shape_function_derivatives_order": 8, 
                            "sbm_parameters": {
                                "is_inner" : false
                            }
                        }
                    ] // element condition list
                }
            }] // iga modeler
            """)
        
        iga_model_part = current_model.CreateModelPart("IgaModelPart")
        iga_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        iga_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)
        
        run_modelers(current_model, modeler_settings)

        # assign the boundary conditions to the trus skin
        skin_model_part = current_model.CreateModelPart("skin_model_part")
        for node in skin_model_part.Nodes:
            node.SetValue(KM.DISPLACEMENT_X, node.X+node.Y)
            node.SetValue(KM.DISPLACEMENT_Y, 123.4)

        self.assertEqual(iga_model_part.GetConditions()[24].Info(), "\"SbmSolidCondition\" #24")
        condition = iga_model_part.GetCondition(24) 

        properties_settings = KM.Parameters("""{
            "properties" : [
                {
                    "model_part_name": "IgaModelPart",
                    "properties_id": 1,
                    "Material": 
                        {
                            "name": "lin_el",
                            "constitutive_law": { "name": "LinearElasticPlaneStrain2DLaw" },
                            "Variables": 
                                        { 
                                            "PENALTY_FACTOR": -1,
                                            "THICKNESS": 1.0,
                                            "YOUNG_MODULUS": 10,
                                            "POISSON_RATIO": 0.3,
                                            "DENSITY": 1
                                        },
                            "Tables": {}
                        }
                }
            ]}
        """)

        KM.ReadMaterialsUtility(current_model).ReadMaterials(properties_settings)

        # Initialize and step the condition
        process_info = iga_model_part.ProcessInfo

        condition.Initialize(process_info)
        
        # Compute local system
        lhs = KM.Matrix()
        rhs = KM.Vector()
        condition.CalculateLocalSystem(lhs, rhs, process_info)

        # Check the RHS
        expected_RHS = [
            74.6524,37.4204,-0.263736,-37.1943,0,0,0,0,19.8499,176.72967,-1.25275,-176.67315,0,0,0,0,-80.0383,77.2433,-0.549451,-77.4882,0,0,0,0,-12.3542,6.16138,-0.043956,-6.19906,0,0,0,0
        ]
        tolerance = 1e-4

        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

        expected_LHS = [
            -0.0762363,-0.0812402,0.131868,0,0,0,0,0,-0.0520261,0.299941,0.626374,0,0,0,0,0,0.107315,0.340757,0.274725,0,0,0,0,0,0.0209478,0.0433673,0.021978,0,0,0,0,0
        ]
        for i in range(rhs.Size()):
            self.assertAlmostEqual(lhs[0, i], expected_LHS[i], delta=tolerance)


    def testSbmLoadSolidCondition(self):
        current_model = KM.Model()
        skin_model_part_outer_initial = current_model.CreateModelPart("skinModelPart_outer_initial")
        skin_model_part_outer_initial.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        skin_model_part_outer_initial.CreateNewProperties(1)
        skin_model_part_outer_initial.CreateNewNode(1, 0.0, 0.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(2, 2.0, 0.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(3, 2.0, 2.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(4, 0.0, 2.0, 0.0)
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_model_part_outer_initial.GetProperties()[1])
        
        skin_model_part = current_model.CreateModelPart("skin_model_part")
        skin_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        
        modeler_settings = KM.Parameters("""
        [
            {
            "modeler_name": "NurbsGeometryModelerSbm",
            "Parameters": {
                    "model_part_name" : "IgaModelPart",
                    "lower_point_xyz": [0.0,0.0,0.0],
                    "upper_point_xyz": [2.0,2.0,0.0],
                    "lower_point_uvw": [0.0,0.0,0.0],
                    "upper_point_uvw": [2.0,2.0,0.0],
                    "polynomial_order" : [3, 3],
                    "number_of_knot_spans" : [5,5],
                    "lambda_outer": 0.5,
                    "number_of_inner_loops": 0,
                    "skin_model_part_outer_initial_name": "skinModelPart_outer_initial",           
                    "skin_model_part_name": "skin_model_part"
                }
            },
            {
                "modeler_name": "IgaModelerSbm",
                "Parameters": {
                    "echo_level": 0,
                    "skin_model_part_name": "skin_model_part",
                    "analysis_model_part_name": "IgaModelPart",
                    "element_condition_list": [
                        {
                            "geometry_type": "SurfaceEdge",
                            "iga_model_part": "SBM_Support_outer",
                            "type": "condition",
                            "name": "SbmLoadSolidCondition",
                            "shape_function_derivatives_order": 4, 
                            "sbm_parameters": {
                                "is_inner" : false
                            }
                        }
                    ] // element condition list
                }
            }] // iga modeler
            """)
        
        iga_model_part = current_model.CreateModelPart("IgaModelPart")
        iga_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        iga_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)
        
        run_modelers(current_model, modeler_settings)

        # assign the boundary conditions to the trus skin
        skin_model_part = current_model.CreateModelPart("skin_model_part")
        for node in skin_model_part.Nodes:
            node.SetValue(KM.FORCE_X, node.X+node.Y)
            node.SetValue(KM.FORCE_Y, 123.4)

        self.assertEqual(iga_model_part.GetConditions()[24].Info(), "\"SbmLoadSolidCondition\" #24")
        condition = iga_model_part.GetCondition(24) 

        properties_settings = KM.Parameters("""{
            "properties" : [
                {
                    "model_part_name": "IgaModelPart",
                    "properties_id": 1,
                    "Material": 
                        {
                            "name": "lin_el",
                            "constitutive_law": { "name": "LinearElasticPlaneStrain2DLaw" },
                            "Variables": 
                                        { 
                                            "PENALTY_FACTOR": -1,
                                            "THICKNESS": 1.0,
                                            "YOUNG_MODULUS": 10,
                                            "POISSON_RATIO": 0.3,
                                            "DENSITY": 1
                                        },
                            "Tables": {}
                        }
                }
            ]}
        """)

        KM.ReadMaterialsUtility(current_model).ReadMaterials(properties_settings)

        # Initialize and step the condition
        process_info = iga_model_part.ProcessInfo

        condition.Initialize(process_info)
        
        # Compute local system
        lhs = KM.Matrix()
        rhs = KM.Vector()
        condition.CalculateLocalSystem(lhs, rhs, process_info)

        # Check the RHS
        expected_RHS = [
            0.00261224,1.2894,0,0,0,0,0,0,0.0124082,6.12467,0,0,0,0,0,0,0.00544218,2.68626,0,0,0,0,0,0,0.000435374,0.214901,0,0,0,0,0,0
        ]
        
        tolerance = 1e-4

        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

        expected_LHS = [
            -0.074175824175824, -0.049450549450549, 0.074175824175824, 0.0, 0.0, 0.0, 0.0, 0.0, -0.055631868131868, 0.058281004709576, 0.055631868131868, 0.0, 0.0, 0.0, 0.0, 0.0, 0.109203296703297, 0.001766091051805, -0.109203296703297, 0.0, 0.0, 0.0, 0.0, 0.0, 0.020604395604396, -0.010596546310832, -0.020604395604396, 0.0, 0.0, 0.0, 0.0, 0.0
        ]

        for i in range(rhs.Size()):
            self.assertAlmostEqual(lhs[0, i], expected_LHS[i], delta=tolerance)

if __name__ == '__main__':
    KratosUnittest.main()