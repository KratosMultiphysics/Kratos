# Kratos Imports
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

from KratosMultiphysics.StructuralMechanicsApplication.RBE3_process import ApplyRbe3Process


class TestRBE3(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KM.ROTATION)
        mp.AddNodalSolutionStepVariable(KM.REACTION)
        mp.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)


    def _add_dofs(self,mp):
        # Adding the dofs AND their corresponding reaction!
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X,mp)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y,mp)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z,mp)
        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X,mp)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y,mp)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z,mp)

    def _create_nodes(self, mp):
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 100.0, 0.0, 0.0)
        mp.CreateNewNode(3, 200.0, 0.0, 0.0)
        mp.CreateNewNode(4, 300.0, 0.0, 0.0)
        mp.CreateNewNode(5, 400.0, 0.0, 0.0)
        mp.CreateNewNode(6, 500.0, 0.0, 0.0)
        mp.CreateNewNode(7, 600.0, 0.0, 0.0)

        mp.CreateNewNode(8, 0.0, 0.0, -100.0)
        mp.CreateNewNode(9, 100.0, 0.0, -100.0)
        mp.CreateNewNode(10, 200.0, 0.0, -100.0)
        mp.CreateNewNode(11, 300.0, 0.0, -100.0)
        mp.CreateNewNode(12, 400.0, 0.0, -100.0)
        mp.CreateNewNode(13, 500.0, 0.0, -100.0)
        mp.CreateNewNode(14, 600.0, 0.0, -100.0)

        mp.CreateNewNode(15, 0.0, 0.0, -200.0)
        mp.CreateNewNode(16, 100.0, 0.0, -200.0)
        mp.CreateNewNode(17, 200.0, 0.0, -200.0)
        mp.CreateNewNode(18, 300.0, 0.0, -200.0)
        mp.CreateNewNode(19, 400.0, 0.0, -200.0)
        mp.CreateNewNode(20, 500.0, 0.0, -200.0)
        mp.CreateNewNode(21, 600.0, 0.0, -200.0)

        mp.CreateNewNode(22, 0.0, 0.0, -300.0)
        mp.CreateNewNode(23, 100.0, 0.0, -300.0)
        mp.CreateNewNode(24, 200.0, 0.0, -300.0)
        mp.CreateNewNode(25, 300.0, 0.0, -300.0)
        mp.CreateNewNode(26, 400.0, 0.0, -300.0)
        mp.CreateNewNode(27, 500.0, 0.0, -300.0)
        mp.CreateNewNode(28, 600.0, 0.0, -300.0)

        mp.CreateNewNode(29, 0.0, 0.0, -400.0)
        mp.CreateNewNode(30, 100.0, 0.0, -400.0)
        mp.CreateNewNode(31, 200.0, 0.0, -400.0)
        mp.CreateNewNode(32, 300.0, 0.0, -400.0)
        mp.CreateNewNode(33, 400.0, 0.0, -400.0)
        mp.CreateNewNode(34, 500.0, 0.0, -400.0)
        mp.CreateNewNode(35, 600.0, 0.0, -400.0)

        mp.CreateNewNode(36, 0.0, 0.0, -500.0)
        mp.CreateNewNode(37, 100.0, 0.0, -500.0)
        mp.CreateNewNode(38, 200.0, 0.0, -500.0)
        mp.CreateNewNode(39, 300.0, 0.0, -500.0)
        mp.CreateNewNode(40, 400.0, 0.0, -500.0)
        mp.CreateNewNode(41, 500.0, 0.0, -500.0)
        mp.CreateNewNode(42, 600.0, 0.0, -500.0)

        mp.CreateNewNode(43, 0.0, 0.0, -600.0)
        mp.CreateNewNode(44, 100.0, 0.0, -600.0)
        mp.CreateNewNode(45, 200.0, 0.0, -600.0)
        mp.CreateNewNode(46, 300.0, 0.0, -600.0)
        mp.CreateNewNode(47, 400.0, 0.0, -600.0)
        mp.CreateNewNode(48, 500.0, 0.0, -600.0)
        mp.CreateNewNode(49, 600.0, 0.0, -600.0)

        mp.CreateNewNode(50, 0.0, 0.0, -700.0)
        mp.CreateNewNode(51, 100.0, 0.0, -700.0)
        mp.CreateNewNode(52, 200.0, 0.0, -700.0)
        mp.CreateNewNode(53, 300.0, 0.0, -700.0)
        mp.CreateNewNode(54, 400.0, 0.0, -700.0)
        mp.CreateNewNode(55, 500.0, 0.0, -700.0)
        mp.CreateNewNode(56, 600.0, 0.0, -700.0)

        mp.CreateNewNode(57, 0.0, 0.0, -800.0)
        mp.CreateNewNode(58, 100.0, 0.0, -800.0)
        mp.CreateNewNode(59, 200.0, 0.0, -800.0)
        mp.CreateNewNode(60, 300.0, 0.0, -800.0)
        mp.CreateNewNode(61, 400.0, 0.0, -800.0)
        mp.CreateNewNode(62, 500.0, 0.0, -800.0)
        mp.CreateNewNode(63, 600.0, 0.0, -800.0)

        mp.CreateNewNode(64, 0.0, 0.0, -900.0)
        mp.CreateNewNode(65, 100.0, 0.0, -900.0)
        mp.CreateNewNode(66, 200.0, 0.0, -900.0)
        mp.CreateNewNode(67, 300.0, 0.0, -900.0)
        mp.CreateNewNode(68, 400.0, 0.0, -900.0)
        mp.CreateNewNode(69, 500.0, 0.0, -900.0)
        mp.CreateNewNode(70, 600.0, 0.0, -900.0)

        mp.CreateNewNode(71, 0.0, 0.0, -1000.0)
        mp.CreateNewNode(72, 100.0, 0.0, -1000.0)
        mp.CreateNewNode(73, 200.0, 0.0, -1000.0)
        mp.CreateNewNode(74, 300.0, 0.0, -1000.0)
        mp.CreateNewNode(75, 400.0, 0.0, -1000.0)
        mp.CreateNewNode(76, 500.0, 0.0, -1000.0)
        mp.CreateNewNode(77, 600.0, 0.0, -1000.0)

        mp.CreateNewNode(1000, 470.0, 150.0, -1120.0)

    def _create_elements(self, mp, element_name):
       mp.CreateNewElement(element_name, 1, [1, 2, 9, 8], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 2, [2, 3, 10, 9], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 3, [3, 4, 11, 10], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 4, [4, 5, 12, 11], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 5, [5, 6, 13, 12], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 6, [6, 7, 14, 13], mp.GetProperties()[1])

       mp.CreateNewElement(element_name, 7, [8, 9, 16, 15], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 8, [9, 10, 17, 16], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 9, [10, 11, 18, 17], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 10, [11, 12, 19, 18], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 11, [12, 13, 20, 19], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 12, [13, 14, 21, 20], mp.GetProperties()[1])

       mp.CreateNewElement(element_name, 13, [15, 16, 23, 22], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 14, [16, 17, 24, 23], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 15, [17, 18, 25, 24], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 16, [18, 19, 26, 25], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 17, [19, 20, 27, 26], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 18, [20, 21, 28, 27], mp.GetProperties()[1])

       mp.CreateNewElement(element_name, 19, [22, 23, 30, 29], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 20, [23, 24, 31, 30], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 21, [24, 25, 32, 31], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 22, [25, 26, 33, 32], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 23, [26, 27, 34, 33], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 24, [27, 28, 35, 34], mp.GetProperties()[1])

       mp.CreateNewElement(element_name, 25, [29, 30, 37, 36], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 26, [30, 31, 38, 37], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 27, [31, 32, 39, 38], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 28, [32, 33, 40, 39], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 29, [33, 34, 41, 40], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 30, [34, 35, 42, 41], mp.GetProperties()[1])

       mp.CreateNewElement(element_name, 31, [36, 37, 44, 43], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 32, [37, 38, 45, 44], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 33, [38, 39, 46, 45], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 34, [39, 40, 47, 46], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 35, [40, 41, 48, 47], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 36, [41, 42, 49, 48], mp.GetProperties()[1])

       mp.CreateNewElement(element_name, 37, [43, 44, 51, 50], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 38, [44, 45, 52, 51], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 39, [45, 46, 53, 52], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 40, [46, 47, 54, 53], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 41, [47, 48, 55, 54], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 42, [48, 49, 56, 55], mp.GetProperties()[1])

       mp.CreateNewElement(element_name, 43, [50, 51, 58, 57], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 44, [51, 52, 59, 58], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 45, [52, 53, 60, 59], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 46, [53, 54, 61, 60], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 47, [54, 55, 62, 61], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 48, [55, 56, 63, 62], mp.GetProperties()[1])

       mp.CreateNewElement(element_name, 49, [57, 58, 65, 64], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 50, [58, 59, 66, 65], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 51, [59, 60, 67, 66], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 52, [60, 61, 68, 67], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 53, [61, 62, 69, 68], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 54, [62, 63, 70, 69], mp.GetProperties()[1])

       mp.CreateNewElement(element_name, 55, [64, 65, 72, 71], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 56, [65, 66, 73, 72], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 57, [66, 67, 74, 73], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 58, [67, 68, 75, 74], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 59, [68, 69, 76, 75], mp.GetProperties()[1])
       mp.CreateNewElement(element_name, 60, [69, 70, 77, 76], mp.GetProperties()[1])


    def _apply_dirichlet_BCs(self,mp):
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_X, True, mp.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Y, True, mp.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Z, True, mp.Nodes)
        KM.VariableUtils().ApplyFixity(KM.ROTATION_X, True, mp.Nodes)
        KM.VariableUtils().ApplyFixity(KM.ROTATION_Y, True, mp.Nodes)
        KM.VariableUtils().ApplyFixity(KM.ROTATION_Z, True, mp.Nodes)

    def _apply_neumann_BCs(self,mp,node, applied_load):
        node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD,0,applied_load) #0 = step
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[node.Id],mp.GetProperties()[1])   #ID

    def _apply_neumann_BCs_M(self,mp, node, applied_moment):
        node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_MOMENT,0,applied_moment)
        mp.CreateNewCondition("PointMomentCondition3D1N",2,[node.Id],mp.GetProperties()[1])   

    def _apply_material_properties(self,mp):
        mp.CreateNewProperties(1)
        mp.GetProperties()[1].SetValue(KM.YOUNG_MODULUS,210E+3)
        mp.GetProperties()[1].SetValue(KM.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KM.THICKNESS,2.0)
        mp.GetProperties()[1].SetValue(KM.DENSITY,7.850E-9)

    def _solve(self,mp):
        linear_solver = KM.SkylineLUFactorizationSolver()
        builder_and_solver = KM.ResidualBasedBlockBuilderAndSolver(
            linear_solver
        )
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()

        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True
        strategy = KM.ResidualBasedLinearStrategy(mp,
                                                  scheme,
                                                  builder_and_solver,
                                                  compute_reactions,
                                                  reform_step_dofs,
                                                  calculate_norm_dx,
                                                  move_mesh_flag)
        strategy.SetEchoLevel(0)
        strategy.Initialize()
        strategy.Check()
        strategy.Solve()

    def _check_slave_disp(self,node,displacement_results, rotation_results):
        disp = node.GetSolutionStepValue(KM.DISPLACEMENT)
        rot = node.GetSolutionStepValue(KM.ROTATION)
        for component in range(3):
            self.assertIsClose(disp[component], 
                               displacement_results[component],
                               abs_tol=0.0001, 
                               rel_tol=0.07, 
                               msg= None 
            )
            self.assertIsClose(
                rot[component],
                rotation_results[component],
                abs_tol=0.0001, 
                rel_tol=0.07, 
                msg= None
            )

    def execute_RBE3_test(self, current_model, element_name, applied_load, applied_moment, displacement_results, rotation_results):
        mp = current_model.CreateModelPart("Structure")
        model = current_model
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)

        if element_name == "MITCThickShellElement3D4N":
            cl = KM.StructuralMechanicsApplication.ReissnerMindlinShellElasticConstitutiveLaw()
            mp.GetProperties()[1].SetValue(KM.CONSTITUTIVE_LAW, cl)

        self._create_nodes(mp)
        self._add_dofs(mp)
        self._create_elements(mp,element_name)

        connected = mp.CreateSubModelPart("RBE3_connected") 
        connected.AddNodes([65,67,69,72,74,76])

        reference = mp.CreateSubModelPart("RBE3_reference")
        reference.AddNodes([1000])

        bcs_dirichlet = mp.CreateSubModelPart("BoundaryCondtionsDirichlet")
        bcs_dirichlet.AddNodes([1,2,3,4,5,6,7])

        self._apply_dirichlet_BCs(bcs_dirichlet)

        rbe3_settings = KM.Parameters(r"""
        {
            "model_part_name": "Structure",
            "connected_sub_model_part": "Structure.RBE3_connected",
            "reference_sub_model_part": "Structure.RBE3_reference",
            "constraint_id_start": 1,
            "constrained_dofs": ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]
        }""")

        rbe3_process = ApplyRbe3Process(
            current_model,
            rbe3_settings
        )

        rbe3_process.ExecuteInitialize()

        ref_node = mp.GetNode(1000)

        self._apply_neumann_BCs(mp, ref_node, applied_load)
        self._apply_neumann_BCs_M(mp, ref_node, applied_moment)
        self._solve(mp)


        for node_id, expected_displacement, expected_rotation in zip(
            [65,67,69,72,74,76], displacement_results, rotation_results
        ):
            self._check_slave_disp(
                mp.GetNode(node_id),
                expected_displacement,
                expected_rotation
            )

    def test_RBE3_panel_fx(self):
        self._test_RBE3_panel(
            applied_load=[1000.0, 0.0, 0.0],
            applied_moment=[0.0, 0.0, 0.0],
            displacement_results=[
                [5.39E-02,  1.74E+02, -1.59E-02],  # Slave 65
                [5.37E-02, -6.50E-10, -1.03E-16],  # Slave 67
                [5.39E-02, -1.74E+02,  1.59E-02],  # Slave 69
                [6.28E-02,  1.96E+02, -1.64E-02],  # Slave 72
                [6.24E-02, -7.50E-10, -1.06E-16],  # Slave 74
                [6.28E-02, -1.96E+02,  1.64E-02],  # Slave 76
            ],
            rotation_results=[
                [ 2.14E-01, -8.50E-05, -8.31E-01],  # Slave 65
                [-9.92E-13, -8.24E-05, -8.91E-01],  # Slave 67
                [-2.14E-01, -8.50E-05, -8.31E-01],  # Slave 69
                [ 2.19E-01, -8.50E-05, -9.05E-01],  # Slave 72
                [-1.00E-12, -8.25E-05, -1.01E+00],  # Slave 74
                [-2.19E-01, -8.50E-05, -9.05E-01],  # Slave 76
            ]
    )
    def test_RBE3_panel_fy(self):
        self._test_RBE3_panel(
            applied_load=[0.0, 1000.0, 0.0],
            applied_moment=[0.0, 0.0, 0.0],
            displacement_results=[
                [-1.12E-17, 3.56E+03, 7.75E-18],  # Slave 65
                [-1.10E-17, 3.79E+03, 6.73E-18],  # Slave 67
                [-1.06E-17, 3.96E+03, 5.53E-18],  # Slave 69
                [-1.19E-17, 4.24E+03, 8.09E-18],  # Slave 72
                [-1.16E-17, 4.49E+03, 7.34E-18],  # Slave 74
                [-1.19E-17, 4.69E+03, 6.46E-18],  # Slave 76
            ],
            rotation_results=[
                [6.69E+00, 0.00E+00, 1.21E+00],   # Slave 65
                [6.92E+00, 0.00E+00, 1.01E+00],   # Slave 67
                [7.17E+00, 1.30E-20, 6.75E-01],   # Slave 69
                [6.85E+00, 0.00E+00, 1.26E+00],   # Slave 72
                [7.09E+00, 0.00E+00, 1.15E+00],   # Slave 74
                [7.35E+00, 1.25E-20, 7.95E-01],   # Slave 76
            ]
        )

    def test_RBE3_panel_fz(self):
        self._test_RBE3_panel(
            applied_load=[0.0, 0.0, 1000.0],
            applied_moment=[0.0, 0.0, 0.0],
            displacement_results=[
                [8.81E-03, 6.93E+02, -5.86E-04],  # Slave 65
                [8.87E-03, 7.03E+02,  3.73E-03],  # Slave 67
                [9.15E-03, 6.93E+02,  7.95E-03],  # Slave 69
                [1.11E-02, 8.56E+02, -6.47E-04],  # Slave 72
                [1.12E-02, 8.65E+02,  4.07E-03],  # Slave 74
                [1.11E-02, 8.56E+02,  8.72E-03],  # Slave 76
            ],
            rotation_results=[
                [1.56E+00, -2.04E-05,  1.02E-01],  # Slave 65
                [1.56E+00, -2.04E-05,  1.04E-12],  # Slave 67
                [1.56E+00, -1.79E-05, -1.02E-01],  # Slave 69
                [1.67E+00, -2.17E-05,  7.95E-02],  # Slave 72
                [1.66E+00, -2.17E-05,  1.07E-12],  # Slave 74
                [1.67E+00, -1.84E-05, -7.95E-02],  # Slave 76
            ]
        )
    def test_RBE3_panel_mx(self):
        self._test_RBE3_panel(
            applied_load=[0.0, 0.0, 0.0],
            applied_moment=[100000.0, 0.0, 0.0],
            displacement_results=[
                [1.63E-18, 4.62E+02, -3.52E-18],  # Slave 65
                [1.48E-18, 4.68E+02, -2.06E-18],  # Slave 67
                [1.44E-18, 4.62E+02, -5.50E-19],  # Slave 69
                [2.21E-18, 5.71E+02, -3.93E-18],  # Slave 72
                [2.27E-18, 5.77E+02, -2.25E-18],  # Slave 74
                [2.20E-18, 5.71E+02, -5.25E-19],  # Slave 76
            ],
            rotation_results=[
                [1.04E+00, 0.00E+00,  6.81E-02],  # Slave 65
                [1.04E+00, 0.00E+00,  6.92E-13],  # Slave 67
                [1.04E+00, 0.00E+00, -6.81E-02],  # Slave 69
                [1.12E+00, 0.00E+00,  5.30E-02],  # Slave 72
                [1.11E+00, 0.00E+00,  7.13E-13],  # Slave 74
                [1.12E+00, 0.00E+00, -5.30E-02],  # Slave 76
            ]
        )
    def test_RBE3_panel_my(self):
        self._test_RBE3_panel(
            applied_load=[0.0, 0.0, 0.0],
            applied_moment=[0.0, 100000.0, 0.0],
            displacement_results=[
                [-5.28E-03, -2.08E-14,  2.51E-03],  # Slave 65
                [-5.22E-03, -1.44E-14,  1.15E-17],  # Slave 67
                [-5.28E-03, -7.18E-15, -2.51E-03],  # Slave 69
                [-6.53E-03, -2.64E-14,  2.75E-03],  # Slave 72
                [-6.57E-03, -1.91E-14,  1.20E-17],  # Slave 74
                [-6.53E-03, -1.11E-14, -2.75E-03],  # Slave 76
            ],
            rotation_results=[
                [-5.20E-17, 1.13E-05, 2.81E-17],  # Slave 65
                [-4.38E-17, 1.20E-05, 3.47E-17],  # Slave 67
                [-3.53E-17, 1.13E-05, 3.67E-17],  # Slave 69
                [-5.80E-17, 1.18E-05, 3.21E-17],  # Slave 72
                [-4.90E-17, 1.27E-05, 3.94E-17],  # Slave 74
                [-4.09E-17, 1.18E-05, 3.85E-17],  # Slave 76
            ]
        ) 

    def test_RBE3_panel_mz(self):
        self._test_RBE3_panel(
            applied_load=[0.0, 0.0, 0.0],
            applied_moment=[0.0, 0.0, 100000.0],
            displacement_results=[
                [-1.85E-18, -1.16E+02,  5.61E-19],  # Slave 65
                [-1.85E-18,  4.33E-10,  0.00E+00],  # Slave 67
                [-1.85E-18,  1.16E+02, -5.61E-19],  # Slave 69
                [-2.19E-18, -1.31E+02,  5.82E-19],  # Slave 72
                [-2.17E-18,  4.99E-10,  0.00E+00],  # Slave 74
                [-2.19E-18,  1.31E+02, -5.82E-19],  # Slave 76
            ],
            rotation_results=[
                [-1.42E-01, 0.00E+00, 5.54E-01],  # Slave 65
                [ 6.61E-13, 0.00E+00, 5.94E-01],  # Slave 67
                [ 1.42E-01, 0.00E+00, 5.54E-01],  # Slave 69
                [-1.46E-01, 0.00E+00, 6.03E-01],  # Slave 72
                [ 6.70E-13, 0.00E+00, 6.74E-01],  # Slave 74
                [ 1.46E-01, 0.00E+00, 6.03E-01],  # Slave 76
            ]
        )


    def _test_RBE3_panel(self, applied_load, applied_moment, displacement_results, rotation_results):
        element_name = "MITCThickShellElement3D4N"
        current_model = KM.Model()
        self.execute_RBE3_test(current_model,
                                element_name,
                                applied_load,
                                applied_moment,
                                displacement_results,
                                rotation_results) 

if __name__ == "__main__":
    KratosUnittest.main()


