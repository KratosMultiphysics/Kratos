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
    def test_contained_refinement_with_boundary_refinement_has_closed_hole(self):
        model = KM.Model()
        body = model.CreateModelPart("IgaModelPart").CreateSubModelPart("Body")
        body.ProcessInfo[KM.DOMAIN_SIZE] = 2
        settings = KM.Parameters('''[
            {
                "modeler_name": "ImportNurbsSbmModeler",
                "Parameters": {
                    "input_filename": "import_nurbs_test/circle_00_22.json",
                    "model_part_name": "outer_initial",
                    "link_layer_to_condition_name": [
                        {"layer_name": "Layer0", "condition_name": "SupportLaplacianCondition"}
                    ]
                }
            },
            {
                "modeler_name": "LocalRefinementModeler",
                "Parameters": {
                    "model_part_name": "IgaModelPart.Body",
                    "refinement_type": "sbm",
                    "base_domain": {
                        "lower_point_uvw": [0.7, 0.7, 0],
                        "upper_point_uvw": [1.3, 1.3, 0]
                    },
                    "refinement_regions": [
                        {
                            "lower_point_uvw": [0.7, 0.94, 0],
                            "upper_point_uvw": [0.94, 1.06, 0],
                            "polynomial_order": [1, 1],
                            "number_of_knot_spans": [16, 16],
                            "lambda_outer": 0.5
                        },
                        {
                            "lower_point_uvw": [1.02, 0.98, 0],
                            "upper_point_uvw": [1.06, 1.02, 0],
                            "polynomial_order": [1, 1],
                            "number_of_knot_spans": [8, 8],
                            "lambda_outer": 0.5
                        }
                    ],
                    "geometry_parameters": {
                        "model_part_name": "IgaModelPart.Body",
                        "skin_model_part_name": "skin_Body",
                        "skin_model_part_outer_initial_name": "outer_initial",
                        "polynomial_order": [1, 1],
                        "number_of_knot_spans": [40, 40],
                        "lambda_inner": 0.0,
                        "lambda_outer": 1.0,
                        "number_of_inner_loops": 0,
                        "number_initial_points_if_importing_nurbs": 100,
                        "gap_approximation_order": 1,
                        "number_internal_divisions": 0,
                        "gap_sbm_type": "interpolation",
                        "gap_element_name": "LaplacianElement",
                        "gap_interface_condition_name": "SupportLaplacianCondition"
                    },
                    "analysis_parameters": {
                        "analysis_model_part_name": "IgaModelPart.Body",
                        "element_condition_list": [
                            {"geometry_type": "GeometrySurface",
                             "iga_model_part": "StructuralAnalysisDomain",
                             "type": "element", "name": "LaplacianElement",
                             "shape_function_derivatives_order": 2}
                        ]
                    }
                }
            }
        ]''')
        run_modelers(model, settings)
        patch = model["IgaModelPart.Body.Patch1"]
        inner = patch.GetSubModelPart("surrogate_inner")
        self.assertGreater(inner.NumberOfConditions(), 0)
        degrees = {}
        incoming = {}
        outgoing = {}
        signed_area = 0.0
        for condition in inner.Conditions:
            for node in condition.GetGeometry():
                degrees[node.Id] = degrees.get(node.Id, 0) + 1
                self.assertLess(node.X, 1.1)
                self.assertGreater(node.X, 0.99)
            # CreateBrepsSbmUtilities reverses entering segments. Check the
            # effective BREP orientation, not the raw scan-line direction.
            a, b = condition.GetGeometry()[0], condition.GetGeometry()[1]
            if condition.Is(KM.BOUNDARY):
                a, b = b, a
            outgoing[a.Id] = outgoing.get(a.Id, 0) + 1
            incoming[b.Id] = incoming.get(b.Id, 0) + 1
            signed_area += 0.5 * (a.X * b.Y - b.X * a.Y)
        self.assertTrue(all(degree == 2 for degree in degrees.values()))
        self.assertEqual(set(incoming), set(degrees))
        self.assertEqual(set(outgoing), set(degrees))
        self.assertTrue(all(degree == 1 for degree in incoming.values()))
        self.assertTrue(all(degree == 1 for degree in outgoing.values()))
        self.assertLess(signed_area, 0.0)  # A hole has clockwise orientation.

        def is_inside(boundary, point):
            crossings = 0
            for condition in boundary.Conditions:
                a, b = condition.GetGeometry()[0], condition.GetGeometry()[1]
                if (a.Y > point.Y) != (b.Y > point.Y):
                    x = a.X + (point.Y - a.Y) * (b.X - a.X) / (b.Y - a.Y)
                    crossings += x > point.X
            return crossings % 2 == 1

        outer = patch.GetSubModelPart("surrogate_outer")
        domain = patch.GetSubModelPart("StructuralAnalysisDomain")
        self.assertGreater(domain.NumberOfElements(), 0)
        for element in domain.Elements:
            center = element.GetGeometry().Center()
            self.assertTrue(is_inside(outer, center), f"Active point outside outer loop: {center}")
            self.assertFalse(is_inside(inner, center), f"Active point inside refinement hole: {center}")

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
        self.assertEqual(comp_mp.NumberOfElements(), 832)
        elem = min(comp_mp.Elements, key=lambda e: e.Id)
        self.assertEqual(elem.Info(), f"LaplacianElement #{elem.Id}")

        # Check first 7 element centers against exact expected values
        expected_elem_centers = [
            (0.890487656721211, 0.7975706522882061, 0.0),
            (0.8964878306162233, 0.8040095715495481, 0.0),
            (0.9043163819433897, 0.8124105630559446, 0.0),
            (0.9103165558384012, 0.8188494823172856, 0.0),
            (0.8900751743628588, 0.8055008015373745, 0.0),
            (0.8945273024261354, 0.8106556074888887, 0.0),
            (0.9003360862494602, 0.8173811897268564, 0.0)
        ]
        for i, expected in enumerate(expected_elem_centers):
            center_point = elems_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)

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
                        "lower_point_uvw": [-1,-1,0.0],
                        "upper_point_uvw": [2,2,0.0],
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
        
        expected_centers = [ (0.15406620770145904, 0.007216934927795088, 0.0), (0.1540237583380274, 0.03550236076110134, 0.0), (0.15396159618270683, 0.07692307692307701, 0.0), (0.15389943402738626, 0.11834379308505268, 0.0), (0.15385698466395462, 0.14662921891835892, 0.0), (0.14662921891835892, 0.38461899488798496, 0.0), (0.11834379308505268, 0.3846331446757955, 0.0), (0.07692307692307701, 0.3846538653942357, 0.0), (0.035502360761101345, 0.3846745861126758, 0.0), (0.007216934927795096, 0.38468873590048636, 0.0) ]

        for i, expected in enumerate(expected_centers):
            center_point = conditions_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)
        
        
        expected_elem_centers = [ (0.14865150109255038, 0.2643431519663906, 0.0), (0.12778533742048212, 0.22906068239034427, 0.0), (0.09527839892947686, 0.17409489728751815, 0.0), (0.05856775491667727, 0.112021092018153, 0.0), (0.026060816425672033, 0.05705530691532692, 0.0), (0.005194652753603742, 0.021772837339280586, 0.0), (0.14865150109255038, 0.281226091570515, 0.0), (0.12778533742048212, 0.25089767618707914, 0.0), (0.09527839892947686, 0.20364970339758534, 0.0), (0.05856775491667726, 0.15029175614329424, 0.0) ]

        for i, expected in enumerate(expected_elem_centers):
            center_point = elems_sorted[i].GetGeometry().Center()
            self.assertAlmostEqual(center_point[0], expected[0], places=12)
            self.assertAlmostEqual(center_point[1], expected[1], places=12)
            self.assertAlmostEqual(center_point[2], expected[2], places=12)


if __name__ == "__main__":
    KratosUnittest.main()
