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
    def test_quadrature_points_gap_sbm_on_circle(self):
        current_model = KM.Model()

        # First import NURBS boundary (outer) using the ImportNurbsSbmModeler
        import_modelers = KM.Parameters(
            """
            [
                {
                    "modeler_name" : "ImportNurbsSbmModeler",
                    "Parameters" : {
                        "input_filename" : "import_nurbs_test/circle_00_22.json",
                        "model_part_name" : "initial_skin_model_part_in",
                        "link_layer_to_condition_name": [
                            {
                                "layer_name" : "Layer0",
                                "condition_name" : "SupportLaplacianCondition" 
                            }
                        ]
                    }
                }
            ]
            """
        )
        
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
                        "echo_level": 0,
                        "lower_point_xyz": [0,0,0.0],
                        "upper_point_xyz": [2,2,0.0],
                        "lower_point_uvw": [0,0,0.0],
                        "upper_point_uvw": [2,2,0.0],
                        "polynomial_order" : [1,1],
                        "number_of_knot_spans" : [63, 63],
                        "number_of_inner_loops": 1,
                        "number_initial_points_if_importing_nurbs": 100,
                        "number_internal_divisions": 0,
                        "gap_approximation_order" : 3,
                        "gap_sbm_type": "default",
                        "skin_model_part_inner_initial_name": "initial_skin_model_part_in",           
                        "skin_model_part_name": "skin_model_part",
                        "gap_element_name": "LaplacianElement",
                        "gap_interface_condition_name": "SupportLaplacianCondition"
                    }
                }
            ]
            """
        )

        run_modelers(current_model, modeler_settings)

        support_mp = current_model.GetModelPart("IgaModelPart.GapInterfaces")
        comp_mp = current_model.GetModelPart("IgaModelPart.GapElements")
        conditions_sorted = sorted(support_mp.Conditions, key=lambda c: c.Id)
        elems_sorted = sorted(comp_mp.Elements, key=lambda e: e.Id)
        self.assertEqual(support_mp.NumberOfConditions(), 260)
        expected_centers = [
            (0.8906146650794393, 0.8321748302233414, 0.0), (0.8902817552078408, 0.8308673230686796, 0.0), (0.8897942474227267, 0.8289526309211348, 0.0), (0.8893067396376129, 0.82703793877359, 0.0), (0.8889738297660142, 0.8257304316189282, 0.0), (0.8900433172033742, 0.7948715461219632, 0.0)
        ]
        for i, expected in enumerate(expected_centers):
            center_point = conditions_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)

        # Check element quadrature point and weight in GAP elements submodel part
        self.assertEqual(comp_mp.NumberOfElements(), 468)
        elem = min(comp_mp.Elements, key=lambda e: e.Id)
        self.assertEqual(elem.Info(), f"LaplacianElement #{elem.Id}")

        # Check first 7 element centers against exact expected values
        expected_elem_centers = [
            (0.8913728268702876, 0.7999212781946125, 0.0), (0.899908859059869, 0.8091745608897168, 0.0), (0.9084448912494507, 0.8184278435848212, 0.0), (0.8903776838425083, 0.8113909863056128, 0.0), (0.8954939153035784, 0.817807523768493, 0.0), (0.900610146764649, 0.8242240612313728, 0.0), (0.8893825408147291, 0.8228606944166132, 0.0)
        ]
        for i, expected in enumerate(expected_elem_centers):
            center_point = elems_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)


    # def test_quadrature_points_and_weights(self):
    #     current_model = KM.Model()

    #     # First import NURBS boundary (outer) using the ImportNurbsSbmModeler
    #     import_modelers = KM.Parameters(
    #         """
    #         [
    #             {
    #                 "modeler_name" : "ImportNurbsSbmModeler",
    #                 "Parameters" : {
    #                     "echo_level" : 0,
    #                     "input_filename" : "import_nurbs_test/square_nurbs.json",
    #                     "model_part_name" : "initial_skin_model_part_out",
    #                     "link_layer_to_condition_name": [
    #                         {"layer_name" : "bottom", "condition_name" : "SbmLaplacianConditionDirichlet"},
    #                         {"layer_name" : "left",   "condition_name" : "SbmLaplacianConditionDirichlet"},
    #                         {"layer_name" : "right",  "condition_name" : "SbmLaplacianConditionNeumann"},
    #                         {"layer_name" : "top",    "condition_name" : "SbmLaplacianConditionNeumann"}
    #                     ]
    #                 }
    #             }
    #         ]
    #         """
    #     )

    #     run_modelers(current_model, import_modelers)

    #     # Target model parts
    #     iga_model_part = current_model.CreateModelPart("IgaModelPart")
    #     iga_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)

    #     skin_model_part = current_model.CreateModelPart("skin_model_part")

    #     modeler_settings = KM.Parameters(
    #         """
    #         [
    #             {
    #                 "modeler_name": "NurbsGeometryModelerGapSbm",
    #                 "Parameters": {
    #                     "echo_level": 0,
    #                     "model_part_name" : "IgaModelPart",
    #                     "lower_point_xyz": [-0.2,-0.2,0.0],
    #                     "upper_point_xyz": [1.0,1.0,0.0],
    #                     "lower_point_uvw": [-0.2,-0.2,0.0],
    #                     "upper_point_uvw": [1.0,1.0,0.0],
    #                     "polynomial_order" : [1, 1],
    #                     "number_of_knot_spans" : [17, 17],
    #                     "number_of_inner_loops": 0,
    #                     "number_initial_points_if_importing_nurbs": 10,
    #                     "number_internal_divisions": 0,
    #                     "skin_model_part_outer_initial_name": "initial_skin_model_part_out",
    #                     "skin_model_part_name": "skin_model_part",
    #                     "gap_element_name": "LaplacianElement",
    #                     "gap_interface_condition_name": "SupportLaplacianCondition",
    #                     "gap_sbm_type": "default"
    #                 }
    #             }
    #         ]
    #         """
    #     )

    #     run_modelers(current_model, modeler_settings)

    #     # Print and check centers of some conditions in GAP interfaces submodel part
    #     support_mp = current_model.GetModelPart("IgaModelPart.GapInterfaces")
    #     self.assertEqual(support_mp.NumberOfConditions(), 168)
    #     conditions_sorted = sorted(support_mp.Conditions, key=lambda c: c.Id)
    #     expected_centers = [
    #         (0.0013259019456383313, 0.0013259019456383313, 0.0),
    #         (0.005882352941176464, 0.005882352941176464, 0.0),
    #         (0.010438803936714596, 0.010438803936714596, 0.0),
    #         (0.010438803936714596, 0.07933283118918327, 0.0),
    #         (0.005882352941176464, 0.06895424836601308, 0.0),
    #         (0.0013259019456383313, 0.05857566554284289, 0.0),
    #     ]
    #     for i, expected in enumerate(expected_centers):
    #         center_point = conditions_sorted[i].GetGeometry().Center()
    #         self.assertAlmostEqual(center_point[0], expected[0], places=12)
    #         self.assertAlmostEqual(center_point[1], expected[1], places=12)
    #         self.assertAlmostEqual(center_point[2], expected[2], places=12)

    #     # Check element quadrature point and weight in GAP elements submodel part
    #     comp_mp = current_model.GetModelPart("IgaModelPart.GapElements")
    #     self.assertEqual(comp_mp.NumberOfElements(), 512)
    #     elem = min(comp_mp.Elements, key=lambda e: e.Id)
    #     self.assertEqual(elem.Info(), f"LaplacianElement #{elem.Id}")

    #     elems_sorted = sorted(comp_mp.Elements, key=lambda e: e.Id)
    #     # Check first 7 element centers against exact expected values
    #     expected_elem_centers = [
    #         (0.010947860656435703, 0.015776462755235827, 0.0),
    #         (0.007882241432852197, 0.012438866491913309, 0.0),
    #         (0.003882464449500949, 0.008084235423886636, 0.0),
    #         (0.000816845225917448, 0.00474663916056412, 0.0),
    #         (0.010947860656435703, 0.03389820105699787, 0.0),
    #         (0.007882241432852197, 0.029539875266290386, 0.0),
    #     ]
    #     for i, expected in enumerate(expected_elem_centers):
    #         center_point = elems_sorted[i].GetGeometry().Center()
    #         self.assertAlmostEqual(center_point[0], expected[0], places=12)
    #         self.assertAlmostEqual(center_point[1], expected[1], places=12)
    #         self.assertAlmostEqual(center_point[2], expected[2], places=12)


if __name__ == "__main__":
    KratosUnittest.main()
