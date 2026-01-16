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
                        "echo_level" : 0,
                        "input_filename" : "import_nurbs_test/square_nurbs.json",
                        "model_part_name" : "initial_skin_model_part_out",
                        "link_layer_to_condition_name": [
                            {"layer_name" : "bottom", "condition_name" : "SbmLaplacianConditionDirichlet"},
                            {"layer_name" : "left",   "condition_name" : "SbmLaplacianConditionDirichlet"},
                            {"layer_name" : "right",  "condition_name" : "SbmLaplacianConditionNeumann"},
                            {"layer_name" : "top",    "condition_name" : "SbmLaplacianConditionNeumann"}
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
                        "model_part_name" : "IgaModelPart",
                        "lower_point_xyz": [-0.2,-0.2,0.0],
                        "upper_point_xyz": [1.0,1.0,0.0],
                        "lower_point_uvw": [-0.2,-0.2,0.0],
                        "upper_point_uvw": [1.0,1.0,0.0],
                        "polynomial_order" : [1, 1],
                        "number_of_knot_spans" : [17, 17],
                        "number_of_inner_loops": 0,
                        "number_initial_points_if_importing_nurbs": 10,
                        "number_internal_divisions": 0,
                        "skin_model_part_outer_initial_name": "initial_skin_model_part_out",
                        "skin_model_part_name": "skin_model_part",
                        "gap_element_name": "LaplacianElement",
                        "gap_interface_condition_name": "SupportLaplacianCondition",
                        "gap_sbm_type": "default"
                    }
                }
            ]
            """
        )

        run_modelers(current_model, modeler_settings)

        # Print and check centers of some conditions in GAP interfaces submodel part
        support_mp = current_model.GetModelPart("IgaModelPart.GapInterfaces")
        self.assertEqual(support_mp.NumberOfConditions(), 168)
        conditions_sorted = sorted(support_mp.Conditions, key=lambda c: c.Id)
        expected_centers = [
            (0.0013259019456383313, 0.0013259019456383313, 0.0),
            (0.005882352941176464, 0.005882352941176464, 0.0),
            (0.010438803936714596, 0.010438803936714596, 0.0),
            (0.010438803936714596, 0.07933283118918327, 0.0),
            (0.005882352941176464, 0.06895424836601308, 0.0),
            (0.0013259019456383313, 0.05857566554284289, 0.0),
        ]
        for i, expected in enumerate(expected_centers):
            center_point = conditions_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)

        # Check element quadrature point and weight in GAP elements submodel part
        comp_mp = current_model.GetModelPart("IgaModelPart.GapElements")
        self.assertEqual(comp_mp.NumberOfElements(), 512)
        elem = min(comp_mp.Elements, key=lambda e: e.Id)
        self.assertEqual(elem.Info(), f"LaplacianElement #{elem.Id}")

        elems_sorted = sorted(comp_mp.Elements, key=lambda e: e.Id)
        # Check first 7 element centers against exact expected values
        expected_elem_centers = [
            (0.010947860656435703, 0.015776462755235827, 0.0),
            (0.007882241432852197, 0.012438866491913309, 0.0),
            (0.003882464449500949, 0.008084235423886636, 0.0),
            (0.000816845225917448, 0.00474663916056412, 0.0),
            (0.010947860656435703, 0.03389820105699787, 0.0),
            (0.007882241432852197, 0.029539875266290386, 0.0),
        ]
        for i, expected in enumerate(expected_elem_centers):
            center_point = elems_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)

    def test_circle_layer_prints(self):
        current_model = KM.Model()

        # Import NURBS boundary (outer) using the ImportNurbsSbmModeler
        import_modelers = KM.Parameters(
            """
            [
                {
                    "modeler_name" : "ImportNurbsSbmModeler",
                    "Parameters" : {
                        "echo_level" : 0,
                        "input_filename" : "import_nurbs_test/circle_nurbs.json",
                        "model_part_name" : "initial_skin_model_part_out",
                        "link_layer_to_condition_name": [
                            {"layer_name" : "Layer0", "condition_name" : "SbmSolidCondition"}
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

        current_model.CreateModelPart("skin_model_part")

        modeler_settings = KM.Parameters(
            """
            [
                {
                    "modeler_name": "NurbsGeometryModelerGapSbm",
                    "Parameters": {
                        "echo_level": 0,
                        "model_part_name" : "IgaModelPart",
                        "lower_point_xyz": [-1.0,-1.0,0.0],
                        "upper_point_xyz": [1.0,1.0,0.0],
                        "lower_point_uvw": [-1.0,-1.0,0.0],
                        "upper_point_uvw": [1.0,1.0,0.0],
                        "polynomial_order" : [1, 1],
                        "number_of_knot_spans" : [17, 17],
                        "number_of_inner_loops": 0,
                        "number_initial_points_if_importing_nurbs": 50,
                        "number_internal_divisions": 0,
                        "skin_model_part_outer_initial_name": "initial_skin_model_part_out",
                        "skin_model_part_name": "skin_model_part",
                        "gap_element_name": "LaplacianElement",
                        "gap_interface_condition_name": "SupportLaplacianCondition",
                        "gap_sbm_type": "default"
                    }
                }
            ]
            """
        )

        run_modelers(current_model, modeler_settings)

        # Print and check centers of some conditions in GAP interfaces submodel part
        support_mp = current_model.GetModelPart("IgaModelPart.GapInterfaces")
        self.assertEqual(support_mp.NumberOfConditions(), 84)
        conditions_sorted = sorted(support_mp.Conditions, key=lambda c: c.Id)
        expected_centers = [
            (-0.17195929666172438, -0.4631787664672619, 0.0),
            (-0.17392843786853043, -0.4407369571414995, 0.0),
            (-0.17589757907533649, -0.4182951478157371, 0.0),
            (-0.18914083728283518, -0.30696294900191184, 0.0),
            (-0.23268204263086079, -0.3511057221678383, 0.0),
            (-0.27622324797888637, -0.3952484953337647, 0.0),
            (0.17589757907533674, -0.4182951478157371, 0.0),
        ]
        for i, expected in enumerate(expected_centers):
            center_point = conditions_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)

        # Check element quadrature point and weight in GAP elements submodel part
        comp_mp = current_model.GetModelPart("IgaModelPart.GapElements")
        self.assertEqual(comp_mp.NumberOfElements(), 448)
        elem = min(comp_mp.Elements, key=lambda e: e.Id)
        self.assertEqual(elem.Info(), f"LaplacianElement #{elem.Id}")

        elems_sorted = sorted(comp_mp.Elements, key=lambda e: e.Id)
        # Check first 7 element centers against exact expected values
        expected_elem_centers = [
            (-0.17668405238067955, -0.40788956359648526, 0.0),
            (-0.17748518308069097, -0.4240023511471369, 0.0),
            (-0.17853043492081228, -0.4450250392510395, 0.0),
            (-0.17933156562082367, -0.4611378268016907, 0.0),
            (-0.17881003813640128, -0.3782471175415114, 0.0),
            (-0.1875899910110614, -0.39816449338046617, 0.0),
            (-0.19904537762686034, -0.4241511067330159, 0.0),
        ]
        for i, expected in enumerate(expected_elem_centers):
            center_point = elems_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)


if __name__ == "__main__":
    KratosUnittest.main()
