import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as DFA
import KratosMultiphysics.IgaApplication as IGA
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
    iga_model_part.SetBufferSize(2)
    iga_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
    iga_model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
    iga_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    skin_model_part_outer_initial = current_model.CreateModelPart("skin_model_part_outer_initial")
    skin_model_part_outer_initial.AddNodalSolutionStepVariable(KM.VELOCITY)
    _create_cube_outer_skin(skin_model_part_outer_initial)

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
                    "number_of_knot_spans" : [4, 4, 4],
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
                            "name": "SbmFluidConditionDirichlet",
                            "shape_function_derivatives_order": 3,
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


def _add_fluid_dofs_and_set_3d_state(model_part):
    for node in model_part.Nodes:
        node.AddDof(KM.VELOCITY_X)
        node.AddDof(KM.VELOCITY_Y)
        node.AddDof(KM.VELOCITY_Z)
        node.AddDof(KM.PRESSURE)

        velocity = node.GetSolutionStepValue(KM.VELOCITY)
        velocity[0] = 0.1 + 0.01 * node.X - 0.02 * node.Y + 0.03 * node.Z
        velocity[1] = -0.05 - 0.02 * node.X + 0.03 * node.Y + 0.01 * node.Z
        velocity[2] = 0.02 + 0.02 * node.X + 0.01 * node.Y - 0.01 * node.Z
        node.SetSolutionStepValue(KM.PRESSURE, 1.0 + 0.2 * node.X - 0.1 * node.Y + 0.05 * node.Z)


class SbmStokes3DTests(KratosUnittest.TestCase):

    def test_SbmFluidConditionDirichlet3D(self):
        current_model, support_model_part = _create_outer_support_model_part()

        self.assertEqual(support_model_part.NumberOfConditions(), 288)

        iga_model_part = current_model.GetModelPart("IgaModelPart")
        skin_model_part = current_model.GetModelPart("skin_model_part")

        _add_fluid_dofs_and_set_3d_state(iga_model_part)

        for node in skin_model_part.Nodes:
            value = KM.Vector(3)
            value[0] = 0.4 + 0.1 * node.X
            value[1] = -0.2 + 0.05 * node.Y
            value[2] = 0.3 - 0.04 * node.Z
            node.SetValue(KM.VELOCITY, value)

        properties_settings = KM.Parameters(
            """
            {
                "properties" : [
                    {
                        "model_part_name": "IgaModelPart.SBM_Support_outer",
                        "properties_id": 1,
                        "Material": {
                            "name": "fluid",
                            "constitutive_law": { "name": "Newtonian3DLaw" },
                            "Variables": {
                                "PENALTY_FACTOR": 100.0,
                                "DENSITY": 1.0,
                                "DYNAMIC_VISCOSITY": 0.8
                            },
                            "Tables": {}
                        }
                    }
                ]
            }
            """
        )
        KM.ReadMaterialsUtility(current_model).ReadMaterials(properties_settings)

        condition = next(iter(support_model_part.Conditions))

        process_info = iga_model_part.ProcessInfo
        condition.Initialize(process_info)

        lhs = KM.Matrix()
        rhs = KM.Vector()
        condition.CalculateLocalSystem(lhs, rhs, process_info)

        self.assertEqual(lhs.Size1(), 32)
        self.assertEqual(lhs.Size2(), 32)
        self.assertEqual(rhs.Size(), 32)

        tolerance = 1e-10
        expected_lhs = [
            1.0755238241099545e-01, -1.3433877766671154e-01, -1.3433877766671154e-01, 9.6723633543579930e-02,
           -7.4855835790106610e-02, 4.2641258225238565e-02, 4.2641258225238565e-02, -4.8361816771789970e-02,
            1.4208273396991739e-01, 3.8804597540563604e-02, -1.2706432724027540e-01, 2.5917019497006095e-02,
           -7.6689640728139720e-02, -2.5974591558571910e-02, 5.6959870831847590e-02, -1.2958509748503048e-02]
        expected_rhs = [
            6.7208353624072850e-02, 4.0204230873218914e-03, 1.0759633060402471e-01, -3.3044199858682774e-02,
            8.1012439291299690e-01, -3.7229962760713550e-01, 5.1891613154282820e-01, 1.6522099929341387e-02,
            3.0476817291480820e+00, -1.4087540909603724e+00, 2.1217250757062933e+00, -8.8541666666666680e-03,
           -1.2977644756851523e+00, 6.0515829548018620e-01, -9.0740420451981340e-01, 4.4270833333333340e-03]

        for i, expected_value in enumerate(expected_lhs):
            self.assertAlmostEqual(lhs[0, i], expected_value, delta=tolerance)

        for i, expected_value in enumerate(expected_rhs):
            self.assertAlmostEqual(rhs[i], expected_value, delta=tolerance)


if __name__ == "__main__":
    KratosUnittest.main()
