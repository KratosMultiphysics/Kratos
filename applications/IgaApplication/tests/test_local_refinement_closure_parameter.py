import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest
from test_nurbs_geometry_modeler_gap_sbm import run_modelers


class TestLocalRefinementClosureParameter(KratosUnittest.TestCase):
    def test_closure_has_unit_speed(self):
        model = KM.Model()
        body = model.CreateModelPart("IgaModelPart").CreateSubModelPart("Body")
        body.ProcessInfo[KM.DOMAIN_SIZE] = 2
        settings = KM.Parameters('''[
            {"modeler_name": "ImportNurbsSbmModeler", "Parameters": {
                "input_filename": "import_nurbs_test/circle_00_22.json",
                "model_part_name": "outer_initial",
                "link_layer_to_condition_name": [
                    {"layer_name": "Layer0", "condition_name": "GapSbmLoadSolidCondition"}
                ]
            }},
            {"modeler_name": "LocalRefinementModeler", "Parameters": {
                "model_part_name": "IgaModelPart.Body",
                "refinement_type": "sbm",
                "coupling_conditions_name": "SolidCouplingCondition",
                "base_domain": {
                    "lower_point_uvw": [0.7, 0.7, 0], "upper_point_uvw": [1.3, 1.3, 0]
                },
                "refinement_regions": [{
                    "lower_point_uvw": [0.7, 0.94, 0], "upper_point_uvw": [0.94, 1.06, 0],
                    "polynomial_order": [1, 1], "number_of_knot_spans": [16, 16],
                    "lambda_outer": 0.5
                }],
                "geometry_parameters": {
                    "model_part_name": "IgaModelPart.Body",
                    "skin_model_part_name": "skin_Body",
                    "skin_model_part_outer_initial_name": "outer_initial",
                    "polynomial_order": [1, 1], "number_of_knot_spans": [20, 20],
                    "lambda_inner": 0.0, "lambda_outer": 1.0, "number_of_inner_loops": 0,
                    "number_initial_points_if_importing_nurbs": 100,
                    "gap_approximation_order": 1, "number_internal_divisions": 0,
                    "gap_sbm_type": "interpolation", "gap_element_name": "GapSbmSolidElement",
                    "gap_interface_condition_name": "GapSbmSolidInterfaceCondition"
                },
                "analysis_parameters": {
                    "analysis_model_part_name": "IgaModelPart.Body",
                    "element_condition_list": [{
                        "geometry_type": "GeometrySurface", "iga_model_part": "StructuralAnalysisDomain",
                        "type": "element", "name": "SolidElement", "shape_function_derivatives_order": 2
                    }]
                }
            }}
        ]''')
        run_modelers(model, settings)
        enhanced = model["IgaModelPart.Body.Patch1.LocalRefinementEnhancedBoundaryConditions"]
        self.assertGreater(enhanced.NumberOfConditions(), 0)
        properties = KM.Properties(1)
        properties.SetValue(KM.CONSTITUTIVE_LAW, SMA.LinearElasticPlaneStress2DLaw())
        properties.SetValue(KM.YOUNG_MODULUS, 1000.0)
        properties.SetValue(KM.POISSON_RATIO, 0.3)
        properties.SetValue(KM.THICKNESS, 1.0)
        for condition in enhanced.Conditions:
            self.assertIn('"GapSbmEnhancedLoadSolidCondition"', condition.Info())
            # With direct quadrature weights, a straight closure must be
            # parameterized by arc length, hence have a unit tangent/normal.
            normal = condition.GetGeometry().Normal(0)
            self.assertAlmostEqual(sum(value * value for value in normal), 1.0, places=10)
            condition.Properties = properties
            condition.Initialize(body.ProcessInfo)
            oriented_normal = condition.GetValue(KM.NORMAL)
            self.assertAlmostEqual(sum(normal[i] * oriented_normal[i] for i in range(3)), 1.0, places=10)


if __name__ == "__main__":
    KratosUnittest.main()
