import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication as SMA


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


class TestGapSbmSolid(KratosUnittest.TestCase):
    def test_import_nurbs_and_check_gap_sbm_solid(self):
        current_model = KM.Model()

        # Import NURBS boundary (outer) using the ImportNurbsSbmModeler
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
                            {"layer_name" : "bottom", "condition_name" : "SbmSolidCondition"},
                            {"layer_name" : "left",   "condition_name" : "SbmSolidCondition"},
                            {"layer_name" : "right",  "condition_name" : "SbmSolidCondition"},
                            {"layer_name" : "top",    "condition_name" : "SbmSolidCondition"}
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
        iga_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        iga_model_part.SetBufferSize(2)
        current_model.CreateModelPart("skin_model_part")

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
                        "gap_element_name": "GapSbmSolidElement",
                        "gap_interface_condition_name": "GapSbmSolidInterfaceCondition",
                        "gap_sbm_type": "default"
                    }
                }
            ]
            """
        )

        run_modelers(current_model, modeler_settings)

        # Check GAP interface condition and print its local system
        support_mp = current_model.GetModelPart("IgaModelPart.GapInterfaces")
        self.assertGreater(support_mp.NumberOfConditions(), 0)
        cond = min(support_mp.Conditions, key=lambda c: c.Id)
        self.assertEqual(cond.Info(), f"\"GapSbmSolidInterfaceCondition\" #{cond.Id}")
        # Reuse element properties for condition
        props_if = iga_model_part.CreateNewProperties(2)
        props_if.SetValue(KM.THICKNESS, 1.0)
        props_if.SetValue(KM.YOUNG_MODULUS, 100.0)
        props_if.SetValue(KM.POISSON_RATIO, 0.3)
        law_if = SMA.LinearElasticPlaneStrain2DLaw()
        props_if.SetValue(KM.CONSTITUTIVE_LAW, law_if)
        cond.Properties = props_if
        cond.Initialize(iga_model_part.ProcessInfo)
        lhs_c = KM.Matrix(); rhs_c = KM.Vector()
        cond.CalculateLocalSystem(lhs_c, rhs_c, iga_model_part.ProcessInfo)
        row0_c = [lhs_c[0, j] for j in range(lhs_c.Size2())]
        rhs_list_c = [rhs_c[j] for j in range(rhs_c.Size())]
        # Exact checks for condition LHS first row and RHS
        expected_row0_c = [
            0.5848400711156697,
            -0.1399578865263838,
            -0.735972953909504,
            0.366657210717135,
            -0.05587498133725014,
            -0.01117499626745038,
            0.2070078641310844,
            -0.2155243279233008,
            -0.09525790145658415,
            0.2378743204582017,
            -0.05587498133725002,
            -0.01117499626745046,
            0.6242229912350037,
            -0.3890072032520359,
            -0.4730901084411695,
            0.1623078790612846,
        ]
        expected_rhs_c = [0.0]*16
        self.assertEqual(len(row0_c), len(expected_row0_c))
        self.assertEqual(len(rhs_list_c), len(expected_rhs_c))
        for i, v in enumerate(expected_row0_c):
            self.assertAlmostEqual(row0_c[i], v, places=12)
        for i, v in enumerate(expected_rhs_c):
            self.assertAlmostEqual(rhs_list_c[i], v, places=12)

        # Check GAP elements submodel part
        comp_mp = current_model.GetModelPart("IgaModelPart.GapElements")
        elem = min(comp_mp.Elements, key=lambda e: e.Id)
        self.assertEqual(elem.Info(), f"GapSbmSolidElement #{elem.Id}")
        # Initialize and compute local system (LHS/RHS)
        props = iga_model_part.CreateNewProperties(1)
        props.SetValue(KM.THICKNESS, 1.0)
        props.SetValue(KM.YOUNG_MODULUS, 100.0)
        props.SetValue(KM.POISSON_RATIO, 0.3)
        law = SMA.LinearElasticPlaneStrain2DLaw()
        props.SetValue(KM.CONSTITUTIVE_LAW, law)
        elem.Properties = props
        # Set body force
        bf = KM.Vector(3)
        bf[0] = 2.5
        bf[1] = 0.0
        bf[2] = 0.0
        elem.SetValue(KM.BODY_FORCE, bf)
        # Add DOFs for DISPLACEMENT on nodes
        for node in iga_model_part.Nodes:
            if not node.HasDofFor(KM.DISPLACEMENT_X):
                node.AddDof(KM.DISPLACEMENT_X)
            if not node.HasDofFor(KM.DISPLACEMENT_Y):
                node.AddDof(KM.DISPLACEMENT_Y)
        elem.Initialize(iga_model_part.ProcessInfo)

        lhs_e = KM.Matrix(); rhs_e = KM.Vector()
        elem.CalculateLocalSystem(lhs_e, rhs_e, iga_model_part.ProcessInfo)

        row0_e = [lhs_e[0, j] for j in range(lhs_e.Size2())]
        rhs_list_e = [rhs_e[j] for j in range(rhs_e.Size())]
        # Exact checks for LHS first row and RHS
        expected_row0 = [
            0.5948474848023295,
            0.005212881722353407,
            -0.592636685250874,
            0.2650719481838142,
            0.03581711831206155,
            -0.00300208217089787,
            -0.03802791786351709,
            -0.2672827477352698,
        ]
        expected_rhs = [
            6.753297290135349e-07,
            0.0,
            5.768374773169533e-05,
            0.0,
            4.06939319423068e-08,
            0.0,
            3.475899850876789e-06,
            0.0,
        ]
        self.assertEqual(len(row0_e), len(expected_row0))
        self.assertEqual(len(rhs_list_e), len(expected_rhs))
        for i, v in enumerate(expected_row0):
            self.assertAlmostEqual(row0_e[i], v, places=12)
        for i, v in enumerate(expected_rhs):
            self.assertAlmostEqual(rhs_list_e[i], v, places=12)


if __name__ == "__main__":
    KratosUnittest.main()
