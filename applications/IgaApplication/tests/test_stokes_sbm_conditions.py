import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.IgaApplication
import KratosMultiphysics.FluidDynamicsApplication


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

class SbmStokesTests(KratosUnittest.TestCase):

    def testSbmStokesCondition(self):
        current_model = KM.Model()
        skin_model_part_outer_initial = current_model.CreateModelPart("skinModelPart_outer_initial")
        skin_model_part_outer_initial.AddNodalSolutionStepVariable(KM.VELOCITY)
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
        skin_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        
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
                            "name": "SbmFluidConditionDirichlet",
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
        iga_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        iga_model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        iga_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)
        
        run_modelers(current_model, modeler_settings)

        # assign the boundary conditions to the trus skin
        skin_model_part = current_model.CreateModelPart("skin_model_part")
        for node in skin_model_part.Nodes:
            node.SetValue(KM.VELOCITY_X, node.X+node.Y)
            node.SetValue(KM.VELOCITY_X, 123.4)

        self.assertEqual(iga_model_part.GetConditions()[24].Info(), "\"SbmFluidConditionDirichlet\" #24")
        condition = iga_model_part.GetCondition(24) 

        properties_settings = KM.Parameters("""{
            "properties" : [
                {
                    "model_part_name": "IgaModelPart",
                    "properties_id": 1,
                    "Material": 
                        {
                            "name": "lin_el",
                            "constitutive_law": { "name": "Newtonian2DLaw" },
                            "Variables": 
                                { 
                                    "PENALTY_FACTOR": 1e6,
                                    "DENSITY"           : 1.0,
                                    "DYNAMIC_VISCOSITY" : 1.0
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

        tolerance = 1e-7

        # Check the RHS
        expected_RHS = [1.223927820016326e+07,-1.289404081632652e+01,1.289404081632652e+00,-1.289404081632652e+01,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,1.263591428793468e+08,-3.223510204081629e+00,6.124669387755096e+00,-6.124669387755096e+01,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,8.405053706666660e+07,1.396854421768706e+01,2.686258503401358e+00,-2.686258503401358e+01,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,9.443879700027201e+06,2.149006802721086e+00,2.149006802721086e-01,-2.149006802721086e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00]
        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)
        expected_LHS = [5.230381479591834e+03,2.132653061224488e-02,-1.306122448979591e-03,1.306122448979590e-02,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,5.399891863520406e+04,8.859693877551011e-02,-6.204081632653055e-03,6.204081632653055e-02,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,3.591857727465984e+04,4.261054421768704e-02,-2.721088435374147e-03,2.721088435374147e-02,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,4.035796079931969e+03,4.200680272108838e-03,-2.176870748299317e-04,2.176870748299317e-03,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00]
        for i in range(rhs.Size()):
            self.assertAlmostEqual(lhs[0, i], expected_LHS[i], delta=tolerance)



if __name__ == '__main__':
    KratosUnittest.main()