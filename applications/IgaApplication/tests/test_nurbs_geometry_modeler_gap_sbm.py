import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
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


class TestNurbsGeometryModelerGapSbm(KratosUnittest.TestCase):
    def test_quadrature_points_and_weights(self):
        current_model = KM.Model()

        # First import NURBS boundary (outer) using the ImportNurbsSbmModeler
        import_modelers = KM.Parameters(
            """
            [
                {
                    "modeler_name" : "ImportNurbsSbmModeler",
                    "Parameters" : {
                        "echo_level" : 4,
                        "input_filename" : "import_nurbs_test/square_nurbs.json",
                        "model_part_name" : "initial_skin_model_part_out",
                        "link_layer_to_condition_name": [
                            {"layer_name" : "left",   "condition_name" : "SbmSolidCondition"},
                            {"layer_name" : "right",  "condition_name" : "SbmSolidCondition"},
                            {"layer_name" : "top",    "condition_name" : "SbmLoadSolidCondition"},
                            {"layer_name" : "bottom", "condition_name" : "SbmLoadSolidCondition"}
                        ]
                    }
                }
            ]
            """
        )

        # import_modelers = KM.Parameters(
        #     """
        #     [
        #         {
        #             "modeler_name" : "ImportNurbsSbmModeler",
        #             "Parameters" : {
        #                 "echo_level" : 4,
        #                 "input_filename" : "import_nurbs_test/circle_nurbs.json",
        #                 "model_part_name" : "initial_skin_model_part_out",
        #                 "link_layer_to_condition_name": [
        #                     {"layer_name" : "Layer0", "condition_name" : "SbmSolidCondition"}
        #                 ]
        #             }
        #         }
        #     ]
        #     """
        # )
        
        run_modelers(current_model, import_modelers)

        # Target model parts
        iga_model_part = current_model.CreateModelPart("IgaModelPart")
        iga_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)

        skin_model_part = current_model.CreateModelPart("skin_model_part")

        modeler_settings = KM.Parameters(
            """
            [
                {
                    "modeler_name": "NurbsGeometryModelerGapSbm",
                    "Parameters": {
                        "echo_level": 4,
                        "model_part_name" : "IgaModelPart",
                        "lower_point_xyz": [-1.0,-1.0,0.0],
                        "upper_point_xyz": [2.0,2.0,0.0],
                        "lower_point_uvw": [-1.0,-1.0,0.0],
                        "upper_point_uvw": [2.0,2.0,0.0],
                        "polynomial_order" : [1, 1],
                        "number_of_knot_spans" : [37, 37],
                        "number_of_inner_loops": 0,
                        "number_initial_points_if_importing_nurbs": 50000,
                        "number_internal_divisions": 0,
                        "skin_model_part_outer_initial_name": "initial_skin_model_part_out",
                        "skin_model_part_name": "skin_model_part",
                        "gap_element_name": "CutSbmSolidElement",
                        "gap_interface_condition_name": "CutSbmSolidInterfaceCondition",
                        "gap_sbm_type": "default"
                    }
                }
            ]
            """
        )

        run_modelers(current_model, modeler_settings)

        # Check condition quadrature point and weight in GAP interfaces submodel part
        support_mp = current_model.GetModelPart("IgaModelPart.GapInterfaces")
        self.assertGreater(support_mp.NumberOfConditions(), 0)
        cond = min(support_mp.Conditions, key=lambda c: c.Id)
        self.assertEqual(cond.Info(), f"\"GapSbmSolidInterfaceCondition\" #{cond.Id}")
        ips_c = cond.GetGeometry().IntegrationPoints()
        self.assertGreater(len(ips_c), 0)
        # Extract first GP location and weight
        gp_c = ips_c[0]
        gp_c_coords = [gp_c.X(), gp_c.Y(), gp_c.Z()]
        gp_c_w = gp_c.Weight()
        print("GapSbmSolidInterfaceCondition first GP (xi,eta,zeta):", gp_c_coords)
        print("GapSbmSolidInterfaceCondition first GP weight:", gp_c_w)

        # Check element quadrature point and weight in GAP elements submodel part
        comp_mp = current_model.GetModelPart("IgaModelPart.GapElements")
        self.assertGreater(comp_mp.NumberOfElements(), 0)
        elem = min(comp_mp.Elements, key=lambda e: e.Id)
        self.assertEqual(elem.Info(), f"CutSbmSolidElement #{elem.Id}")
        ips_e = elem.GetGeometry().IntegrationPoints()
        self.assertGreater(len(ips_e), 0)
        gp_e = ips_e[0]
        gp_e_coords = [gp_e.X(), gp_e.Y(), gp_e.Z()]
        gp_e_w = gp_e.Weight()
        print("CutSbmSolidElement first GP (xi,eta,zeta):", gp_e_coords)
        print("CutSbmSolidElement first GP weight:", gp_e_w)


if __name__ == "__main__":
    KratosUnittest.main()
