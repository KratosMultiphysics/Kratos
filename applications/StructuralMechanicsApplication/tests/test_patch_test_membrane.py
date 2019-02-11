from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestPatchTestMembrane(KratosUnittest.TestCase):

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

    def _add_dofs(self,mp):
        # Adding the dofs AND their corresponding reaction!
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

    def _create_nodes_3d3n(self,mp):
        mp.CreateNewNode(1,   0.0000000000,   0.0000000000,   1.0000000000)
        mp.CreateNewNode(2,   0.1666202260,  -0.0055553047,   0.8333333333)
        mp.CreateNewNode(3,   0.0000000000,   0.0000000000,   0.6666666667)
        mp.CreateNewNode(4,   0.3332937753,  -0.0088892271,   1.0000000000)
        mp.CreateNewNode(5,   0.3332937816,  -0.0088887566,   0.6666666667)
        mp.CreateNewNode(6,   0.1666202247,  -0.0055553047,   0.5000000000)
        mp.CreateNewNode(7,   0.5000001491,  -0.0100000000,   0.8333333333)
        mp.CreateNewNode(8,   0.0000000000,   0.0000000000,   0.3333333333)
        mp.CreateNewNode(9,   0.6667065229,  -0.0088892231,   1.0000000000)
        mp.CreateNewNode(10,   0.5000000000,  -0.0100000000,   0.5000000000)
        mp.CreateNewNode(11,   0.3332937816,  -0.0088887566,   0.3333333333)
        mp.CreateNewNode(12,   0.6667065166,  -0.0088887526,   0.6666666667)
        mp.CreateNewNode(13,   0.1666202260,  -0.0055553047,   0.1666666667)
        mp.CreateNewNode(14,   0.8333799231,  -0.0055553007,   0.8333333333)
        mp.CreateNewNode(15,   0.6667065166,  -0.0088887526,   0.3333333333)
        mp.CreateNewNode(16,   0.5000001491,  -0.0100000000,   0.1666666667)
        mp.CreateNewNode(17,   0.8333799243,  -0.0055553007,   0.5000000000)
        mp.CreateNewNode(18,   0.0000000000,   0.0000000000,   0.0000000000)
        mp.CreateNewNode(19,   1.0000000000,   0.0000000000,   1.0000000000)
        mp.CreateNewNode(20,   1.0000000000,   0.0000000000,   0.6666666667)
        mp.CreateNewNode(21,   0.3332937753,  -0.0088892271,   0.0000000000)
        mp.CreateNewNode(22,   0.8333799231,  -0.0055553007,   0.1666666667)
        mp.CreateNewNode(23,   1.0000000000,   0.0000000000,   0.3333333333)
        mp.CreateNewNode(24,   0.6667065229,  -0.0088892231,   0.0000000000)
        mp.CreateNewNode(25,   1.0000000000,   0.0000000000,   0.0000000000)

    def _create_nodes_3d4n(self,mp):
        mp.CreateNewNode(1,  0.0000000000,  0.0000000000,  1.0000000000)
        mp.CreateNewNode(2,  0.0000000000,  0.0000000000,  0.7500000000)
        mp.CreateNewNode(3,  0.2499498474, -0.0075002175,  1.0000000000)
        mp.CreateNewNode(4,  0.2499498571, -0.0074997471,  0.7500000000)
        mp.CreateNewNode(5,  0.0000000000,  0.0000000000,  0.5000000000)
        mp.CreateNewNode(6,  0.5000001491, -0.0100004706,  1.0000000000)
        mp.CreateNewNode(7,  0.2499498571, -0.0074997471,  0.5000000000)
        mp.CreateNewNode(8,  0.5000001491, -0.0100000000,  0.7500000000)
        mp.CreateNewNode(9,  0.5000001491, -0.0100000000,  0.5000000000)
        mp.CreateNewNode(10, 0.0000000000,  0.0000000000,  0.2500000000)
        mp.CreateNewNode(11, 0.7500504508, -0.0075002115,  1.0000000000)
        mp.CreateNewNode(12, 0.2499498568, -0.0074997471,  0.2500000000)
        mp.CreateNewNode(13, 0.7500504414, -0.0074997411,  0.7500000000)
        mp.CreateNewNode(14, 0.5000001491, -0.0100000000,  0.2500000000)
        mp.CreateNewNode(15, 0.7500504414, -0.0074997411,  0.5000000000)
        mp.CreateNewNode(16, 0.0000000000,  0.0000000000,  0.0000000000)
        mp.CreateNewNode(17, 1.0000000000,  0.0000000000,  1.0000000000)
        mp.CreateNewNode(18, 1.0000000000,  0.0000000000,  0.7500000000)
        mp.CreateNewNode(19, 0.2499498474, -0.0075002175,  0.0000000000)
        mp.CreateNewNode(20, 0.7500504414, -0.0074997411,  0.2500000000)
        mp.CreateNewNode(21, 1.0000000000,  0.0000000000,  0.5000000000)
        mp.CreateNewNode(22, 0.5000001491, -0.0100004706,  0.0000000000)
        mp.CreateNewNode(23, 1.0000000000,  0.0000000000,  0.2500000000)
        mp.CreateNewNode(24, 0.7500504508, -0.0075002115,  0.0000000000)
        mp.CreateNewNode(25, 1.0000000000,  0.0000000000,  0.0000000000)

    def _create_elements_3d3n(self,mp):
        element_name = "PreStressMembraneElement3D3N"
        mp.CreateNewElement(element_name, 1, [21, 13, 18], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [11, 13, 21], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 3, [8, 13, 11], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [18, 13,  8], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 5, [24, 16, 21], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 6, [15, 16, 24], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 7, [11, 16, 15], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 8, [21, 16, 11], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 9, [25, 22, 24], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 10, [23, 22, 25], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 11, [15, 22, 23], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 12, [24, 22, 15], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 13, [11, 6, 8], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 14, [5, 6, 11], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 15, [3, 6, 5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 16, [8, 6, 3], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 17, [15, 10, 11], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 18, [12, 10, 15], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 19, [5, 10, 12], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 20, [11,10, 5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 21, [23, 17, 15], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 22, [20, 17, 23], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 23, [12, 17, 20], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 24, [15, 17, 12], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 25, [5, 2, 3], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 26, [4, 2, 5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 27, [1, 2, 4], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 28, [3, 2, 1], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 29, [12, 7, 5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 30, [9, 7, 12], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 31, [4, 7, 9], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 32, [5, 7, 4], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 33, [20, 14, 12], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 34, [19, 14, 20], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 35, [9, 14, 19], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 36, [12, 14, 9], mp.GetProperties()[1])

    def _create_elements_3d4n(self,mp):
        element_name = "PreStressMembraneElement3D4N"
        mp.CreateNewElement(element_name,  1 , [19, 12, 10, 16], mp.GetProperties()[1])
        mp.CreateNewElement(element_name,  2 , [22, 14, 12, 19], mp.GetProperties()[1])
        mp.CreateNewElement(element_name,  3 , [24, 20, 14, 22], mp.GetProperties()[1])
        mp.CreateNewElement(element_name,  4 , [25, 23, 20, 24], mp.GetProperties()[1])
        mp.CreateNewElement(element_name,  5 , [12,  7,  5, 10], mp.GetProperties()[1])
        mp.CreateNewElement(element_name,  6 , [14,  9,  7, 12], mp.GetProperties()[1])
        mp.CreateNewElement(element_name,  7 , [20, 15,  9, 14], mp.GetProperties()[1])
        mp.CreateNewElement(element_name,  8 , [23, 21, 15, 20], mp.GetProperties()[1])
        mp.CreateNewElement(element_name,  9 , [ 7,  4,  2,  5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 10 , [ 9,  8,  4,  7], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 11 , [15, 13,  8,  9], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 12 , [21, 18, 13, 15], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 13 , [ 4,  3,  1,  2], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 14 , [ 8,  6,  3,  4], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 15 , [13, 11,  6,  8], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 16 , [18, 17, 11, 13], mp.GetProperties()[1])

    def _apply_dirichlet_BCs(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)

    def _apply_self_weight(self, mp):
        for node in mp.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_Y, -9.81)
            node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_Z, 0.0)

    def _apply_material_properties(self,mp):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,1000.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.20)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,0.001)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,700.0)
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.RAYLEIGH_ALPHA,0.03)
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.RAYLEIGH_BETA,0.02)

        constitutive_law = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()

        local_axis_1 = KratosMultiphysics.Vector(3)
        local_axis_1[0] = 1.0
        local_axis_1[1] = 0.0
        local_axis_1[2] = 0.0

        local_axis_2= KratosMultiphysics.Vector(3)
        local_axis_2[0] = 0.0
        local_axis_2[1] = 0.0
        local_axis_2[2] = 1.0

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,constitutive_law)
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.PROJECTION_TYPE_COMBO,"planar")
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.PRESTRESS_AXIS_1_GLOBAL,local_axis_1)
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.PRESTRESS_AXIS_2_GLOBAL,local_axis_2)

        prestress = KratosMultiphysics.Vector(3)
        prestress[0]=1e4        #1e4
        prestress[1]=0.0
        prestress[2]=0.0
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.PRESTRESS_VECTOR,prestress)

    def _solve_static(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = StructuralMechanicsApplication.DisplacementAndOtherDoFCriteria(1e-15,1e-15)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 1000
        compute_reactions = False
        reform_step_dofs = False
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp,
                                                                        scheme,
                                                                        linear_solver,
                                                                        convergence_criterion,
                                                                        builder_and_solver,
                                                                        max_iters,
                                                                        compute_reactions,
                                                                        reform_step_dofs,
                                                                        move_mesh_flag)
        strategy.SetEchoLevel(0)
        strategy.Check()
        strategy.Solve()

    def _solve_dynamic(self,mp):

        #define a minimal newton raphson dynamic solver
        damp_factor_m = -0.30
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m)
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-6,1e-9)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 500
        compute_reactions = True
        reform_step_dofs = False
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp,
                                                                        scheme,
                                                                        linear_solver,
                                                                        convergence_criterion,
                                                                        builder_and_solver,
                                                                        max_iters,
                                                                        compute_reactions,
                                                                        reform_step_dofs,
                                                                        move_mesh_flag)
        strategy.SetEchoLevel(0)

        strategy.Check()
        strategy.Solve()


    def _check_static_results(self,node,displacement_results):
        #check that the results are exact on the node
        displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        self.assertAlmostEqual(displacement[0], displacement_results[0], 4)
        self.assertAlmostEqual(displacement[1], displacement_results[1], 4)
        self.assertAlmostEqual(displacement[2], displacement_results[2], 4)

    def _check_dynamic_results(self,node,step,displacement_results):
        displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        self.assertAlmostEqual(displacement, displacement_results[step], 4)

    def _set_and_fill_buffer(self,mp,buffer_size,delta_time):
        # Set buffer size
        mp.SetBufferSize(buffer_size)

        # Fill buffer
        time = mp.ProcessInfo[KratosMultiphysics.TIME]
        time = time - delta_time * (buffer_size)
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for size in range(0, buffer_size):
            step = size - (buffer_size -1)
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            mp.CloneTimeStep(time)

        mp.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False


    def _set_up_system_3d3n(self,current_model):
        mp = current_model.CreateModelPart("Structure")
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_nodes_3d3n(mp)
        self._add_dofs(mp)
        self._create_elements_3d3n(mp)
        self._apply_self_weight(mp)

        #create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet = mp.CreateSubModelPart("BoundaryCondtionsDirichlet")
        bcs_dirichlet.AddNodes([1,3,8,18,19,20,23,25])
        self._apply_dirichlet_BCs(bcs_dirichlet)

        return mp

    def _set_up_system_3d4n(self,current_model):
        mp = current_model.CreateModelPart("Structure")
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_nodes_3d4n(mp)
        self._add_dofs(mp)
        self._create_elements_3d4n(mp)
        self._apply_self_weight(mp)

        #create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet = mp.CreateSubModelPart("BoundaryCondtionsDirichlet")
        bcs_dirichlet.AddNodes([1,2,5,10,16,17,18,21,23,25])
        self._apply_dirichlet_BCs(bcs_dirichlet)

        return mp


    def test_membrane_3d3n_static(self):
        displacement_results = [-4.628753e-12 , -0.04937043 , -6.483677e-12]

        current_model = KratosMultiphysics.Model()

        mp = self._set_up_system_3d3n(current_model)

        self._solve_static(mp)

        self._check_static_results(mp.Nodes[10],displacement_results)

        #self.__post_process(mp)

    def test_membrane_3d4n_static(self):
        displacement_results = [1.73962e-07 , -0.05202096 , -1.243591e-08]

        current_model = KratosMultiphysics.Model()

        mp = self._set_up_system_3d4n(current_model)

        self._solve_static(mp)

        self._check_static_results(mp.Nodes[9],displacement_results)

        #self.__post_process(mp)


    def test_membrane_3d3n_dynamic(self):

        displacement_results = [-0.004145456940147508, -0.016577305470262357, -0.036009405854487204,
        -0.05724634378436276, -0.07431594485756908, -0.08286293670959116,
        -0.08158307966494918, -0.07223491977349739, -0.058530143571647514,
        -0.044671016341598674, -0.03417034853894327]


        current_model = KratosMultiphysics.Model()
        mp = self._set_up_system_3d3n(current_model)

        #time integration parameters
        dt = 0.05
        time = 0.0
        end_time = 0.5
        step = 0

        self._set_and_fill_buffer(mp,2,dt)

        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)

            self._solve_dynamic(mp)
            self._check_dynamic_results(mp.Nodes[10],step-1,displacement_results)

        #self.__post_process(mp)

    def test_membrane_3d4n_dynamic(self):

        displacement_results = [-0.004416597413161373, -0.017672715828946108, -0.0383282878649957,
        -0.060720299929014065, -0.07850778062395564, -0.08727738281567025,
        -0.08586880143129107, -0.07618411850603453, -0.06201094423792496,
        -0.04762424029880159, -0.03666158466786426]


        current_model = KratosMultiphysics.Model()
        mp = self._set_up_system_3d4n(current_model)

        #time integration parameters
        dt = 0.05
        time = 0.0
        end_time = 0.5
        step = 0

        self._set_and_fill_buffer(mp,2,dt)

        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)

            self._solve_dynamic(mp)
            self._check_dynamic_results(mp.Nodes[9],step-1,displacement_results)

        #self.__post_process(mp)


    def __post_process(self, main_model_part):
        from gid_output_process import GiDOutputProcess
        self.gid_output = GiDOutputProcess(main_model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results"       : ["DISPLACEMENT"],
                                                "gauss_point_results" : []
                                            }
                                        }
                                        """)
                                    )

        self.gid_output.ExecuteInitialize()
        self.gid_output.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteInitializeSolutionStep()
        self.gid_output.PrintOutput()
        self.gid_output.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteFinalize()

if __name__ == '__main__':
    KratosUnittest.main()
