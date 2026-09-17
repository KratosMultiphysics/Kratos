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
        expected_RHS = [4.8957074118530586e+7,-1.2894040816326518e+1,1.2894040816326517e+0,-1.2894040816326518e+1,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,5.0543638777730566e+8,-3.2235102040816295e+0,6.1246693877550955e+0,-6.1246693877550960e+1,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,3.3620206767891127e+8,1.3968544217687061e+1,2.6862585034013584e+0,-2.6862585034013584e+1,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,3.7775512353088401e+7,2.1490068027210860e+0,2.1490068027210862e-1,-2.1490068027210860e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0]
        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)
        expected_LHS = [2.0921548571428564e+4,2.1326530612244876e-2,-1.3061224489795905e-3,1.3061224489795905e-2,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,2.1599568999999989e+5,8.8596938775510115e-2,-6.2040816326530551e-3,6.2040816326530551e-2,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,1.4367427721088426e+5,4.2610544217687035e-2,-2.7210884353741473e-3,2.7210884353741471e-2,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,1.6143178095238081e+4,4.2006802721088381e-3,-2.1768707482993174e-4,2.1768707482993175e-3,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0,0.0000000000000000e+0]
        for i in range(rhs.Size()):
            self.assertAlmostEqual(lhs[0, i], expected_LHS[i], delta=tolerance)



if __name__ == '__main__':
    KratosUnittest.main()
