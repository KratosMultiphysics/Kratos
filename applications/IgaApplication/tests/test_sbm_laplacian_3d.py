import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.ConvectionDiffusionApplication as CDA
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.modeler_factory import KratosModelerFactory


def _run_modelers(current_model, modelers_list):
    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()


def _create_cube_outer_skin(model_part, min_coord=0.25, max_coord=1.75):
    model_part.CreateNewProperties(1)
    prop = model_part.GetProperties()[1]

    coords = {
        1: (min_coord, min_coord, min_coord),
        2: (max_coord, min_coord, min_coord),
        3: (max_coord, max_coord, min_coord),
        4: (min_coord, max_coord, min_coord),
        5: (min_coord, min_coord, max_coord),
        6: (max_coord, min_coord, max_coord),
        7: (max_coord, max_coord, max_coord),
        8: (min_coord, max_coord, max_coord),
    }

    for node_id, (x, y, z) in coords.items():
        model_part.CreateNewNode(node_id, x, y, z)

    triangles = {
        1: [1, 3, 2],
        2: [1, 4, 3],
        3: [5, 6, 7],
        4: [5, 7, 8],
        5: [1, 2, 6],
        6: [1, 6, 5],
        7: [4, 8, 7],
        8: [4, 7, 3],
        9: [1, 5, 8],
        10: [1, 8, 4],
        11: [2, 7, 6],
        12: [2, 3, 7],
    }

    for cond_id, node_ids in triangles.items():
        model_part.CreateNewCondition("SurfaceCondition3D3N", cond_id, node_ids, prop)


def _create_outer_support_model_part():
    current_model = KM.Model()
    iga_model_part = current_model.CreateModelPart("IgaModelPart")
    iga_model_part.CreateNewProperties(0)
    iga_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    skin_model_part = current_model.CreateModelPart("skin_model_part_outer_initial")
    _create_cube_outer_skin(skin_model_part)

    modeler_settings = KM.Parameters(
        """
        [
            {
                "modeler_name": "NurbsGeometryModelerSbm",
                "Parameters": {
                    "model_part_name" : "IgaModelPart",
                    "lower_point_xyz": [0.0, 0.0, 0.0],
                    "upper_point_xyz": [2.0, 2.0, 2.0],
                    "lower_point_uvw": [0.0, 0.0, 0.0],
                    "upper_point_uvw": [2.0, 2.0, 2.0],
                    "polynomial_order" : [1, 1, 1],
                    "number_of_knot_spans" : [2, 2, 2],
                    "lambda_outer": 0.5,
                    "number_of_inner_loops": 0,
                    "skin_model_part_outer_initial_name": "skin_model_part_outer_initial",
                    "skin_model_part_name": "skin_model_part",
                    "echo_level": 0
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
                            "name": "SbmLaplacianConditionDirichlet",
                            "shape_function_derivatives_order": 2,
                            "sbm_parameters": {
                                "is_inner" : false
                            }
                        }
                    ]
                }
            }
        ]
        """
    )

    _run_modelers(current_model, modeler_settings)
    return current_model, current_model.GetModelPart("IgaModelPart.SBM_Support_outer")


class TestSbmLaplacian3D(KratosUnittest.TestCase):

    def test_laplacian_element_3d_local_system(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")

        model_part.AddNodalSolutionStepVariable(KM.CONDUCTIVITY)
        model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(KM.HEAT_FLUX)
        model_part.SetBufferSize(2)
        model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

        settings = KM.ConvectionDiffusionSettings()
        settings.SetDiffusionVariable(KM.CONDUCTIVITY)
        settings.SetUnknownVariable(KM.TEMPERATURE)
        settings.SetVolumeSourceVariable(KM.HEAT_FLUX)
        model_part.ProcessInfo.SetValue(KM.CONVECTION_DIFFUSION_SETTINGS, settings)

        nodes = KM.NodesVector()
        coordinates = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 1.0),
            (0.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
        ]

        for node_id, (x, y, z) in enumerate(coordinates, start=1):
            node = model_part.CreateNewNode(node_id, x, y, z)
            node.AddDof(KM.TEMPERATURE)
            nodes.append(node)

        knot_u = KM.Vector(4)
        knot_v = KM.Vector(4)
        knot_w = KM.Vector(4)
        for knot_vector in (knot_u, knot_v, knot_w):
            knot_vector[0] = 0.0
            knot_vector[1] = 0.0
            knot_vector[2] = 1.0
            knot_vector[3] = 1.0

        volume = KM.NurbsVolumeGeometry(nodes, 1, 1, 1, knot_u, knot_v, knot_w)
        quadrature_point_geometries = KM.GeometriesVector()
        volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2)

        prop = model_part.CreateNewProperties(0)
        prop.SetValue(KM.CONDUCTIVITY, 1.0)

        element = model_part.CreateNewElement("LaplacianElement", 1, quadrature_point_geometries[0], prop)
        element.SetValue(KM.HEAT_FLUX, 1.0)

        lhs = KM.Matrix()
        rhs = KM.Vector()
        element.Initialize(model_part.ProcessInfo)
        element.CalculateLocalSystem(lhs, rhs, model_part.ProcessInfo)

        expected_lhs = [
            1.4508545031536990e-01,
            -2.2444797274783875e-02,
            -2.2444797274783875e-02,
            -2.2444797274783868e-02,
            -2.2444797274783868e-02,
            -2.2444797274783868e-02,
            -2.2444797274783868e-02,
            -1.0416666666666668e-02,
        ]
        expected_rhs = [
            6.1320326520293006e-02,
            1.6430731970725269e-02,
            1.6430731970725269e-02,
            4.4026013626080654e-03,
            1.6430731970725269e-02,
            4.4026013626080663e-03,
            4.4026013626080663e-03,
            1.1796734797069916e-03,
        ]

        tolerance = 1.0e-12
        self.assertEqual(lhs.Size1(), len(expected_lhs))
        self.assertEqual(rhs.Size(), len(expected_rhs))

        for i, expected_value in enumerate(expected_lhs):
            self.assertAlmostEqual(lhs[0, i], expected_value, delta=tolerance)

        for i, expected_value in enumerate(expected_rhs):
            self.assertAlmostEqual(rhs[i], expected_value, delta=tolerance)

    def test_sbm_support_outer_3d_normals(self):
        current_model, support_model_part = _create_outer_support_model_part()
        self.assertTrue(current_model.HasModelPart("IgaModelPart.SBM_Support_outer"))

        self.assertEqual(support_model_part.NumberOfConditions(), 48)

        tolerance = 1.0e-12
        for condition in support_model_part.Conditions:
            center = condition.GetGeometry().Center()
            normal = condition.GetGeometry().Calculate(KM.NORMAL)

            if abs(center[0]) < tolerance:
                expected_normal = (-1.0, 0.0, 0.0)
            elif abs(center[0] - 2.0) < tolerance:
                expected_normal = (1.0, 0.0, 0.0)
            elif abs(center[1]) < tolerance:
                expected_normal = (0.0, -1.0, 0.0)
            elif abs(center[1] - 2.0) < tolerance:
                expected_normal = (0.0, 1.0, 0.0)
            elif abs(center[2]) < tolerance:
                expected_normal = (0.0, 0.0, -1.0)
            else:
                expected_normal = (0.0, 0.0, 1.0)

            self.assertAlmostEqual(normal[0], expected_normal[0], delta=tolerance)
            self.assertAlmostEqual(normal[1], expected_normal[1], delta=tolerance)
            self.assertAlmostEqual(normal[2], expected_normal[2], delta=tolerance)


if __name__ == "__main__":
    KratosUnittest.main()
