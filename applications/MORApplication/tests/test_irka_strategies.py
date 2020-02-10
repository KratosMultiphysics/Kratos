from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.MORApplication as MOR
import KratosMultiphysics.EigenSolversApplication as EigenSolversApplication

from KratosMultiphysics.python_linear_solver_factory import ConstructSolver

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
eigen_solvers_application_available = kratos_utils.CheckIfApplicationsAvailable("EigenSolversApplication")

class IRKATests(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        mp.AddNodalSolutionStepVariable(MOR.COMPONENT_OUTPUT)

    # def _solve_eigen(self,mp,echo=0):
    #     feast_system_solver_settings = KratosMultiphysics.Parameters("""{ }""")

    #     eigensolver_settings = KratosMultiphysics.Parameters("""
    #     {
    #             "perform_stochastic_estimate": false,
    #             "lambda_min": 1.0e-2,
    #             "lambda_max": 25.0,
    #             "search_dimension": 7
    #     }
    #     """)
    #     if not (hasattr(KratosMultiphysics.ExternalSolversApplication,"PastixComplexSolver")):
    #         self.skipTest('"PastixComplexSolver" not available')

    #     feast_system_solver = ExternalSolversApplication.PastixComplexSolver(feast_system_solver_settings)
    #     eigen_solver = ExternalSolversApplication.FEASTSolver(eigensolver_settings, feast_system_solver)
    #     builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(eigen_solver)

    #     eigen_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
    #     eig_strategy = StructuralMechanicsApplication.EigensolverStrategy(mp,
    #                                                                 eigen_scheme,
    #                                                                 builder_and_solver)

    #     eig_strategy.SetEchoLevel(echo)
    #     eig_strategy.Solve()

    # def _setup_harmonic_solver(self,mp,echo=0):
    #     builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(KratosMultiphysics.LinearSolver())
    #     eigen_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
    #     harmonic_strategy = StructuralMechanicsApplication.HarmonicAnalysisStrategy(mp, eigen_scheme, builder_and_solver, False)
    #     harmonic_strategy.SetEchoLevel(echo)

    #     return harmonic_strategy


    def _add_dofs(self,mp):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,mp)

    @KratosUnittest.skipUnless(eigen_solvers_application_available,"Missing required application: EigenSolversApplication")
    def test_complex_irka(self):
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart('mp')
        self._add_variables(mp)

        input_file = 'cantilever_plate'
        material_settings = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                "materials_filename" : "StructuralMaterials_damped.json"
            }
        }
        """)

        KratosMultiphysics.ModelPartIO(input_file).ReadModelPart(mp)
        KratosMultiphysics.ReadMaterialsUtility(material_settings, model)
        self._add_dofs(mp)
        # print(Model)

        mp.SetBufferSize(2)
        mp.Check(mp.ProcessInfo)


        # add load
        for node in model.GetModelPart('mp.PointLoad3D_tipright').Nodes:
            node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD,0,[0,0,-1000])

        for node in model.GetModelPart('mp.DISPLACEMENT_left').Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
            node.Fix(KratosMultiphysics.ROTATION_X)
            node.Fix(KratosMultiphysics.ROTATION_Y)
            node.Fix(KratosMultiphysics.ROTATION_Z)

        output_process_settings = KratosMultiphysics.Parameters("""
            {
                "build_output_structure": true,
                "output_structure_type": "vector",
                "model_part_names" : ["GENERIC_right"],
                "output_variable_names" : ["DISPLACEMENT_Z"]
            }
            """)

        solver_type = 'skyline_lu_complex'
        # solver_type = 'pardiso_ldlt_complex'
        cplx_linear_solver_settings = KratosMultiphysics.Parameters('{ "solver_type" : "' + solver_type + '" }')
        cplx_linear_solver = ConstructSolver(cplx_linear_solver_settings)
        print(cplx_linear_solver)

        solver_type = 'complex_dense_col_piv_householder_qr'
        #solver_type = 'complex_dense_colpivhouseholderqr'
        settings = KratosMultiphysics.Parameters('{ "solver_type" : "EigenSolversApplication.' + solver_type + '" }')
        #cplx_linear_solver = ConstructSolver(settings)
        cplx_dense_linear_solver = EigenSolversApplication.ComplexDenseLinearSolverFactory().Create(settings)

        # dummy_linear_solver = KratosMultiphysics.LinearSolverFactory().Create(dummy_linear_solver_settings)
        #builder_and_solver = MOR.SystemMatrixBuilderAndSolver(dummy_linear_solver,output_process_settings)
        #scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        scheme = MOR.MatrixBuilderScheme()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(KratosMultiphysics.LinearSolver())
        move_mesh_flag = False
        off_strategy = MOR.MorSecondOrderComplexIrkaStrategy(mp,
                                                        scheme,
                                                        builder_and_solver,
                                                        cplx_linear_solver,
                                                        [15j, 25j, 45j],#, 75, 90],
                                                        20,
                                                        1e-12,
                                                        move_mesh_flag)

        off_strategy.SetEchoLevel(3)

        #print(linear_solver)

        strategy = MOR.MorComplexOnlineStrategy(mp,
                                                    #scheme,
                                                    cplx_dense_linear_solver,
                                                    #builder_and_solver,
                                                    off_strategy,
                                                    False)
                                                    #[8, 50, 1500, 3000],
                                                    #move_mesh_flag)

        freq = 10
        max_freq = 30
        df = 10
        step = 1

        # prepare output
        result = []
        n_nodes = len(model.GetModelPart('mp.GENERIC_right').Nodes)
        for cond in model.GetModelPart('mp.GENERIC_right').Conditions:
            cond.SetValue(MOR.COMPONENT_OUTPUT,[0,0,1])

        while freq <= max_freq:
            print('Frequency step: ' + str(freq))

            mp.CloneTimeStep(freq)
            mp.ProcessInfo[MOR.FREQUENCY] = freq
            # vtk_output.ExecuteInitializeSolutionStep()
            # p.ExecuteInitializeSolutionStep()

            strategy.Solve()
            # off_strategy.EchoInfo()
            disp_z = 0
            disp_z = strategy.GetScalarResult()
            # for node in Model.GetModelPart('mp.GENERIC_right').Nodes:
            #     disp_z = disp_z + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0)
            # for node in Model.GetModelPart('mp.plate_top').Nodes:
            #     disp_z = disp_z + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0)
            # for node in Model.GetModelPart('mp.abh_top').Nodes:
            #     disp_z = disp_z + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0)

            # result.update( {freq : disp_z/n_nodes} )
            result.append([freq, disp_z/n_nodes])
            # result = np.vstack([result, np.asarray([freq, disp_z/n_nodes])])

            # result = np.vstack([result, np.asarray([freq, strategy.GetScalarResult()])])
            # print(n_nodes)
            freq = freq + df

            # if step%20 == 0:
            #     vtk_output.PrintOutput()
            # vtk_output.ExecuteFinalizeSolutionStep()
            # p.ExecuteFinalizeSolutionStep()

            step = step + 1

        print(result)
        import csv
        # with open('result.csv', 'w', newline='') as csvfile:
        #     writer = csv.writer(csvfile)
        #     writer.writerows(result)
        with open('result.csv', newline='') as csvfile:
            expected_result = csv.reader(csvfile)
            for r, xr in zip(result, expected_result):
                self.assertAlmostEqual(r[1], complex(xr[1]))

        # print(expected_result)



#     @KratosUnittest.skipUnless(external_solvers_application_available,"Missing required application: ExternalSolversApplication")
#     def test_damped_mdof_harmonic(self):
#         import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
#         current_model = KratosMultiphysics.Model()

#         #analytic solution taken from Humar - Dynamics of Structures p. 677

#         #material properties
#         stiffness = 10.0
#         mass = 2.0
#         damping = 0.1

#         #create the model
#         mp = self._create_2dof_geometry(current_model,stiffness, mass, damping)

#         #solve the eigenproblem
#         self._solve_eigen(mp)

#         #perform the harmonic analysis
#         harmonic_solver = self._setup_harmonic_solver(mp)

#         exfreq = 1.0
#         max_exfreq = 2.0
#         df = 0.05

#         while(exfreq <= max_exfreq):
#             mp.CloneTimeStep(exfreq)
#             harmonic_solver.Solve()

#             exfreq_d = exfreq / sqrt(stiffness/mass)

#             disp_x1_expected_complex = 1 / (3*stiffness * complex(0.5 - exfreq_d**2, sqrt(2) * damping * exfreq_d)) \
#                 + 2 / (3*stiffness * complex(2 - exfreq_d**2, 2 * sqrt(2) * damping * exfreq_d))
#             if disp_x1_expected_complex.real < 0:
#                 disp_x1_expected = - abs(disp_x1_expected_complex)
#             else:
#                 disp_x1_expected = abs(disp_x1_expected_complex)

#             disp_x2_expected_complex = 2 / (3*stiffness * complex(0.5 - exfreq_d**2, sqrt(2) * damping * exfreq_d)) \
#                 - 2 / (3*stiffness * complex(2 - exfreq_d**2, 2 * sqrt(2) * damping * exfreq_d))
#             if disp_x2_expected_complex.real < 0:
#                 disp_x2_expected = - abs(disp_x2_expected_complex)
#             else:
#                 disp_x2_expected = abs(disp_x2_expected_complex)

#             #test displacement
#             self.assertAlmostEqual(mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0), \
#                 disp_x1_expected,delta=1e-5)
#             self.assertAlmostEqual(mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0), \
#                 disp_x2_expected,delta=1e-5)

#             #test phase angle
#             self.assertAlmostEqual(mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.REACTION_X,0), \
#                 phase(disp_x1_expected_complex),delta=1e-5)
#             self.assertAlmostEqual(mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.REACTION_X,0), \
#                 phase(disp_x2_expected_complex),delta=1e-5)

#             exfreq = exfreq + df

# class HarmonicAnalysisTestsWithHDF5(KratosUnittest.TestCase):
#     @KratosUnittest.skipUnless(hdf5_application_available and eigen_solvers_application_available,"Missing required application: HDF5Application, EigenSolversApplication")
#     def test_harmonic_mdpa_input(self):

#         with KratosUnittest.WorkFolderScope(".",__file__):
#             #run simulation and write to hdf5 file
#             model = KratosMultiphysics.Model()
#             project_parameter_file_name = "harmonic_analysis_test/harmonic_analysis_test_eigenproblem_parameters.json"
#             with open(project_parameter_file_name,'r') as parameter_file:
#                 project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
#             test = structural_mechanics_analysis.StructuralMechanicsAnalysis(model,project_parameters)
#             test.Run()

#             #start new simulation and read from hdf5 file
#             model = KratosMultiphysics.Model()
#             project_parameter_file_name = "harmonic_analysis_test/harmonic_analysis_test_parameters.json"
#             with open(project_parameter_file_name,'r') as parameter_file:
#                 project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
#             test = structural_mechanics_analysis.StructuralMechanicsAnalysis(model,project_parameters)
#             test.Run()

#             # remove hdf5 file
#             kratos_utils.DeleteFileIfExisting("harmonic_analysis_test/eigen_results.h5")
#             kratos_utils.DeleteFileIfExisting("harmonic_analysis_test/harmonic_analysis_test.time")

if __name__ == '__main__':
    KratosUnittest.main()
