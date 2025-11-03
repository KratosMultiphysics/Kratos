import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest


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


class GapSbmModelerTests(KratosUnittest.TestCase):
    def test_gap_sbm_interface_condition_print(self):
        if not KM.kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
            self.skipTest("StructuralMechanicsApplication not available")

        current_model = KM.Model()

        # Create initial outer skin on a square [-0.5,0.5]^2
        skin_outer_init = current_model.CreateModelPart("initial_skin_model_part_out")
        skin_outer_init.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        skin_outer_init.CreateNewProperties(1)
        skin_outer_init.CreateNewNode(1, -0.5, -0.5, 0.0)
        skin_outer_init.CreateNewNode(2,  0.5, -0.5, 0.0)
        skin_outer_init.CreateNewNode(3,  0.5,  0.5, 0.0)
        skin_outer_init.CreateNewNode(4, -0.5,  0.5, 0.0)
        skin_outer_init.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_outer_init.GetProperties()[1])
        skin_outer_init.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_outer_init.GetProperties()[1])
        skin_outer_init.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_outer_init.GetProperties()[1])
        skin_outer_init.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_outer_init.GetProperties()[1])

        # Target model part
        iga_model_part = current_model.CreateModelPart("IgaModelPart")
        iga_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        iga_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)

        modeler_settings = KM.Parameters(
            """
            [
                {
                    "modeler_name": "NurbsGeometryModelerSbm",
                    "Parameters": {
                        "echo_level": 4,
                        "model_part_name" : "IgaModelPart",
                        "lower_point_xyz": [-1,-1,0.0],
                        "upper_point_xyz": [ 1, 1,0.0],
                        "lower_point_uvw": [-1,-1,0.0],
                        "upper_point_uvw": [ 1, 1,0.0],
                        "polynomial_order" : [1, 1],
                        "number_of_knot_spans" : [7, 7],
                        "number_of_inner_loops": 0,
                        "skin_model_part_outer_initial_name": "initial_skin_model_part_out",
                        "skin_model_part_name": "skin_model_part"
                    }
                },
                {
                    "modeler_name": "IgaModelerSbm",
                    "Parameters": {
                        "echo_level": 4,
                        "skin_model_part_name": "skin_model_part",
                        "analysis_model_part_name": "IgaModelPart",
                        "element_condition_list": [
                            {
                                "geometry_type": "GeometrySurface",
                                "iga_model_part": "ComputationalDomain",
                                "type": "element",
                                "name": "GapSbmSolidElement",
                                "shape_function_derivatives_order": 3,
                                "variables": [
                                    { "variable_name": "BODY_FORCE", "value": ["30.0","0.0","0.0"] }
                                ]
                            },
                            {
                                "geometry_type": "SurfaceEdge",
                                "iga_model_part": "SBM_Support_outer",
                                "type": "condition",
                                "name": "GapSbmSolidInterfaceCondition",
                                "shape_function_derivatives_order": 3,
                                "sbm_parameters": { "is_inner": false }
                            }
                        ]
                    }
                }
            ]
            """
        )

        run_modelers(current_model, modeler_settings)

        # Set material properties
        properties_settings = KM.Parameters(
            """{
                "properties" : [
                    {
                        "model_part_name": "IgaModelPart",
                        "properties_id": 1,
                        "Material": {
                            "name": "lin_el",
                            "constitutive_law": { "name": "LinearElasticPlaneStrain2DLaw" },
                            "Variables": { 
                                "THICKNESS": 1.0,
                                "YOUNG_MODULUS": 100.0,
                                "POISSON_RATIO": 0.3,
                                "DENSITY": 1.0
                            },
                            "Tables": {}
                        }
                    }
                ]}
            """
        )
        KM.ReadMaterialsUtility(current_model).ReadMaterials(properties_settings)

        process_info = iga_model_part.ProcessInfo

        # Print one interface condition LHS/RHS for copy-paste
        support_mp = current_model.GetModelPart("IgaModelPart.SBM_Support_outer")

        print(support_mp.NumberOfConditions())
        
        # check number of conditions SBM_Support_outer
        self.assertEqual(support_mp.NumberOfConditions(), 36)
        # Take the first quadrature point (smallest Id)
        cond = min(support_mp.Conditions, key=lambda c: c.Id)


        print(cond)
        print(cond.Info())
        # Verify correct condition name via Info()
        self.assertEqual(cond.Info(), f"\"GapSbmSolidInterfaceCondition\" #{cond.Id}")

        # Initialize and compute
        cond.Initialize(process_info)

        print("end initialize")

        lhs = KM.Matrix(); rhs = KM.Vector()
        
        cond.CalculateLocalSystem(lhs, rhs, process_info)

        print(lhs)
        print(rhs)

        # Pretty-print first row and RHS for baseline capture
        row0 = [lhs[0, j] for j in range(lhs.Size2())]
        rhs_list = [rhs[j] for j in range(rhs.Size())]
        print("GapSbmSolidInterfaceCondition LHS row 0:")
        print(row0)
        print("GapSbmSolidInterfaceCondition RHS:")
        print(rhs_list)

        # Note: GapSbmSolidElement requires NEIGHBOUR_GEOMETRIES, which is not exposed in Python.
        # We therefore avoid assembling it here. The element LHS/RHS is covered in C++ tests.

if __name__ == "__main__":
    KratosUnittest.main()
