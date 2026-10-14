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
    # Checks Gap-SBM quadrature results for a circular NURBS boundary.
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

        modeler_settings = KM.Parameters(
            """
            [
                {
                    "modeler_name": "NurbsGeometryModelerGapSbm",
                    "Parameters": {
                        "echo_level": 0,
                        "lower_point_xyz": [0,0,0.0],
                        "upper_point_xyz": [2,2,0.0],
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
        self.assertEqual(support_mp.NumberOfConditions(), 156)

        expected_centers = [
            (0.890495535127507, 0.8317069460367353, 0.0),
            (0.8897942474227267, 0.8289526309211348, 0.0),
            (0.8890929597179466, 0.8261983158055343, 0.0),
            (0.8916624079670692, 0.7965836568056914, 0.0),
            (0.9011935831844301, 0.8066624166158509, 0.0),
            (0.910724758401791, 0.8167411764260103, 0.0)
        ]
        for i, expected in enumerate(expected_centers):
            center_point = conditions_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)

        # Check element quadrature point and weight in GAP elements submodel part
        # 52 gap patches with a 2x2 Gauss rule per patch.
        self.assertEqual(comp_mp.NumberOfElements(), 208)
        elem = min(comp_mp.Elements, key=lambda e: e.Id)
        self.assertEqual(elem.Info(), f"LaplacianElement #{elem.Id}")

        # Check first 7 element centers against exact expected values
        expected_elem_centers = [
            (0.893071317016737, 0.805014318931252, 0.0),
            (0.9044979231610226, 0.8177314673762837, 0.0),
            (0.8902896860333689, 0.8210355193076179, 0.0),
            (0.8941167350031858, 0.8274487894588172, 0.0),
            (1.109710313966631, 0.8210355193076179, 0.0),
            (1.1058832649968149, 0.8274487894588176, 0.0),
            (1.1069286829832627, 0.805014318931252, 0.0)
        ]
        for i, expected in enumerate(expected_elem_centers):
            center_point = elems_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)

    # Checks Gap-SBM quadrature results for a square boundary with layers.
    def test_quadrature_points_gap_sbm_on_square_layers(self):
        current_model = KM.Model()

        import_modelers = KM.Parameters(
            """
            [
                {
                    "modeler_name" : "ImportNurbsSbmModeler",
                    "Parameters" : {
                        "input_filename" : "import_nurbs_test/square_nurbs.json",
                        "model_part_name" : "initial_skin_model_part_out",
                        "link_layer_to_condition_name": [
                            {
                                "layer_name" : "top",
                                "condition_name" : "SupportLaplacianCondition"
                            },
                            {
                                "layer_name" : "left",
                                "condition_name" : "SupportLaplacianCondition"
                            },
                            {
                                "layer_name" : "bottom",
                                "condition_name" : "SupportLaplacianCondition"
                            },
                            {
                                "layer_name" : "right",
                                "condition_name" : "SupportLaplacianCondition"
                            }
                        ]
                    }
                }
            ]
            """
        )
        
        run_modelers(current_model, import_modelers)
        
        iga_model_part = current_model.CreateModelPart("IgaModelPart")
        iga_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)
        current_model.CreateModelPart("skin_model_part")

        modeler_settings = KM.Parameters(
            """
            [
                {
                    "modeler_name": "NurbsGeometryModelerGapSbm",
                    "Parameters": {
                        "echo_level": 0,
                        "model_part_name" : "IgaModelPart",
                        "lower_point_xyz": [-1,-1,0.0],
                        "upper_point_xyz": [2,2,0.0],
                        "polynomial_order" : [2,2],
                        "number_of_knot_spans" : [13, 13],
                        "number_of_inner_loops": 0,
                        "number_initial_points_if_importing_nurbs": 2000,
                        "number_internal_divisions": 0,
                        "gap_approximation_order" : 1,
                        "gap_sbm_type": "interpolation",
                        "skin_model_part_outer_initial_name": "initial_skin_model_part_out",
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
        
        # The equidistant lower-left corner projects onto the left skin edge.
        expected_centers = [
            (0.007216934927795088, 0.154066207701459, 0.0),
            (0.03550236076110134, 0.15402375833802737, 0.0),
            (0.07692307692307701, 0.1539615961827068, 0.0),
            (0.11834379308505268, 0.15389943402738623, 0.0),
            (0.14662921891835892, 0.1538569846639546, 0.0),
            (0.14662921891835892, 0.38461899488798496, 0.0),
            (0.11834379308505268, 0.3846331446757955, 0.0),
            (0.07692307692307701, 0.3846538653942357, 0.0),
            (0.035502360761101345, 0.3846745861126758, 0.0),
            (0.007216934927795096, 0.38468873590048636, 0.0)
        ]

        for i, expected in enumerate(expected_centers):
            center_point = conditions_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)
        
        
        expected_elem_centers = [
            (0.13650743609549892, 0.1798782964794003, 0.0),
            (0.07692307692307707, 0.1799609991129664, 0.0),
            (0.017338717750655198, 0.18004370174653236, 0.0),
            (0.1365074360954989, 0.26924811662221576, 0.0),
            (0.07692307692307707, 0.26930773078847126, 0.0),
            (0.0173387177506552, 0.26936734495472664, 0.0),
            (0.13650743609549892, 0.3586179367650313, 0.0),
            (0.07692307692307707, 0.35865446246397636, 0.0),
            (0.017338717750655205, 0.35869098816292105, 0.0),
            (0.8634925639045016, 0.3586179367650314, 0.0)
        ]

        for i, expected in enumerate(expected_elem_centers):
            center_point = elems_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)


if __name__ == "__main__":
    KratosUnittest.main()
