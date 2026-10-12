import KratosMultiphysics
import KratosMultiphysics.IgaApplication
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

class TestImportNurbsModeler(KratosUnittest.TestCase):
    def testRectangleCurve2D(self):
        current_model = KratosMultiphysics.Model()
        
        modelers_list = KratosMultiphysics.Parameters(
        """ [{
        "modeler_name" : "ImportNurbsSbmModeler",
        "Parameters" : {
            "input_filename" : "import_nurbs_test/square_nurbs.json",
            "model_part_name" : "square_nurbs_model_part",
            "link_layer_to_condition_name": [
            {
                "layer_name" : "left",
                "condition_name" : "SbmSolidCondition"
            },
            {
                "layer_name" : "right",
                "condition_name" : "SbmSolidCondition"
            },
            {
                "layer_name" : "top",
                "condition_name" : "SbmLoadSolidCondition"
            },
            {
                "layer_name" : "bottom",
                "condition_name" : "SbmLoadSolidCondition"
            }
            ]
            }
        }] """ )

        run_modelers(current_model, modelers_list)

        nurbs_model_part = current_model.GetModelPart("square_nurbs_model_part")

        nurbs_curve_center = [
            [0.5, 0.0],
            [1.0, 0.5],
            [0.5, 1.0],
            [0.0, 0.5]
        ]
        nurbs_curve_weigth = [
            [1, 1],
            [1, 1],
            [1, 1],
            [1, 1]
        ]

        i = 0
        for geom in nurbs_model_part.Geometries:
            self.assertEqual(geom.Weights().Size(), 2)
            
            for i_cp in range(0, 2):
                self.assertEqual(geom.Weights()[i_cp], nurbs_curve_weigth[i][i_cp])
            self.assertEqual(geom.Center().X, nurbs_curve_center[i][0])
            self.assertEqual(geom.Center().Y, nurbs_curve_center[i][1])
            i += 1
    

    def testComplexCurve2D(self):
        current_model = KratosMultiphysics.Model()
        
        modelers_list = KratosMultiphysics.Parameters(
        """ [{
        "modeler_name" : "ImportNurbsSbmModeler",
        "Parameters" : {
            "input_filename" : "import_nurbs_test/complex_2d_nurbs.json",
            "model_part_name" : "square_nurbs_model_part",
            "link_layer_to_condition_name": [
            {
                "layer_name" : "Layer0",
                "condition_name" : "SbmLoadSolidCondition"
            },
            {
                "layer_name" : "Layer1",
                "condition_name" : "SbmSolidCondition"
            }
            ]
            }
        }] """ )

        run_modelers(current_model, modelers_list)

        nurbs_model_part = current_model.GetModelPart("square_nurbs_model_part")

        nurbs_curve_center = [
            [-2.416340349142, 4.774201445766],
            [-9.314935000000, 1.458561000000],
            [-3.778479772647, -1.133425102702],
            [3.116020000000, -1.149169000000]
        ]
        nurbs_curve_weigth = [
            [1,1,1,1,1,1,1],
            [1,1],
            [1,1,1,1,1,1],
            [1, 1]
        ]
        nurbs_curve_control_points = [7,2,6,2]

        i = 0
        for geom in nurbs_model_part.Geometries:
            self.assertEqual(geom.Weights().Size(), nurbs_curve_control_points[i])
            
            for i_cp in range(0, 2):
                self.assertEqual(geom.Weights()[i_cp], nurbs_curve_weigth[i][i_cp])
            self.assertAlmostEqual(geom.Center().X, nurbs_curve_center[i][0])
            self.assertAlmostEqual(geom.Center().Y, nurbs_curve_center[i][1])
            i += 1


    def _make_square_modeler_settings(self):
        return KratosMultiphysics.Parameters(
            """ [{
                "modeler_name": "ImportNurbsSbmModeler",
                "Parameters": {
                    "input_filename": "import_nurbs_test/square_nurbs.json",
                    "model_part_name": "square_nurbs_model_part",
                    "link_layer_to_condition_name": [
                        {"layer_name": "bottom", "condition_name": "SbmSolidCondition"},
                        {"layer_name": "right", "condition_name": "SbmSolidCondition"},
                        {"layer_name": "top", "condition_name": "SbmSolidCondition"},
                        {"layer_name": "left", "condition_name": "SbmSolidCondition"}
                    ]
                }
            }] """)

    def test_missing_layer_mapping_for_nurbs_curve(self):
        modelers_list = self._make_square_modeler_settings()
        layers = modelers_list[0]["Parameters"]["link_layer_to_condition_name"]
        layers[0]["layer_name"].SetString("other")

        message = (
            'NURBS geometry layer "bottom" has no entry in '
            '"link_layer_to_condition_name"'
        )
        with self.assertRaisesRegex(RuntimeError, message):
            run_modelers(KratosMultiphysics.Model(), modelers_list)

    def test_missing_condition_name_for_nurbs_layer(self):
        modelers_list = self._make_square_modeler_settings()
        layers = modelers_list[0]["Parameters"]["link_layer_to_condition_name"]
        layers[0].RemoveValue("condition_name")

        with self.assertRaisesRegex(RuntimeError, 'layer "bottom" has no "condition_name"'):
            run_modelers(KratosMultiphysics.Model(), modelers_list)

    def test_empty_condition_name_for_nurbs_layer(self):
        modelers_list = self._make_square_modeler_settings()
        layers = modelers_list[0]["Parameters"]["link_layer_to_condition_name"]
        layers[0]["condition_name"].SetString("")

        message = 'Empty "condition_name" for NURBS geometry layer "bottom"'
        with self.assertRaisesRegex(RuntimeError, message):
            run_modelers(KratosMultiphysics.Model(), modelers_list)

    def test_imports_curve_with_unregistered_condition_name(self):
        modelers_list = self._make_square_modeler_settings()
        layers = modelers_list[0]["Parameters"]["link_layer_to_condition_name"]
        layers[0]["condition_name"].SetString("UnknownSbmCondition")

        current_model = KratosMultiphysics.Model()
        run_modelers(current_model, modelers_list)

        self.assertEqual(
            current_model.GetModelPart("square_nurbs_model_part").NumberOfGeometries(), 4)


if __name__ == '__main__':
    KratosUnittest.main()
