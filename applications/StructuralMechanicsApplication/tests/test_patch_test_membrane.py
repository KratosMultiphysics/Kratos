from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

class BasePatchTestMembrane(KratosUnittest.TestCase):

    def _add_variables(self,mp,explicit_dynamics=False):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        if explicit_dynamics:
            mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_VELOCITY)
            mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.FRACTIONAL_ACCELERATION)
            mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.FRACTIONAL_ANGULAR_ACCELERATION)
            mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
            mp.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
            mp.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)
            mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_ANGULAR_VELOCITY)
            mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.NODAL_INERTIA)
            mp.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENT_RESIDUAL)

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
        element_name = "MembraneElement3D3N"
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
        element_name = "MembraneElement3D4N"
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

    def _apply_dirichlet_BCs(self,mp,fix_type='all'):
        if fix_type=='all':
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)
        elif fix_type=='YZ':
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)
        else:
            print('fix_type: ', fix_type,' not implemented')

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

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,constitutive_law)

        prestress = KratosMultiphysics.Vector(3)
        prestress[0]=1e4        #1e4
        prestress[1]=0.0
        prestress[2]=0.0
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.PRESTRESS_VECTOR,prestress)

    def _solve_static(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-6,1e-6)
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

    def _set_up_system_3d3n(self,current_model,explicit_dynamics=False):
        mp = current_model.CreateModelPart("Structure")
        mp.SetBufferSize(2)
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        self._add_variables(mp,explicit_dynamics)
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
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

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


class StaticPatchTestMembrane(BasePatchTestMembrane):

    def test_membrane_3d3n_static(self):
        displacement_results = [0.0 ,-0.14387, 0.0]

        current_model = KratosMultiphysics.Model()

        mp = self._set_up_system_3d3n(current_model)

        self._solve_static(mp)

        self._check_static_results(mp.Nodes[10],displacement_results)



    def test_membrane_3d4n_static(self):
        displacement_results = [0.0 , -0.594047 , 0.0]

        current_model = KratosMultiphysics.Model()

        mp = self._set_up_system_3d4n(current_model)

        self._solve_static(mp)

        self._check_static_results(mp.Nodes[9],displacement_results)

        #self.__post_process(mp)

    def test_membrane_wrinkling_law(self):

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("Structure")
        mp.SetBufferSize(2)
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        self._add_variables(mp)

        # add properties and subproperties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,206900000000.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.30)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,0.0001)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,7850.0)
        constitutive_law = StructuralMechanicsApplication.WrinklingLinear2DLaw()
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,constitutive_law)
        sub_constitutive_law = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
        mp.GetProperties()[2].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,sub_constitutive_law)
        mp.GetProperties()[1].AddSubProperties(mp.GetProperties()[2])


        # create nodes
        mp.CreateNewNode(1,   0.0000000000,   1.0000000000,   0.0000000000)
        mp.CreateNewNode(2,   0.0000000000,   0.0000000000,   0.0000000000)
        mp.CreateNewNode(3,   1.0000000000,   1.0000000000,   0.0000000000)
        mp.CreateNewNode(4,   1.0000000000,   0.0000000000,   0.0000000000)

        # add dofs
        self._add_dofs(mp)

        # create element
        element_name = "MembraneElement3D4N"
        mp.CreateNewElement(element_name, 1, [4, 3, 1, 2], mp.GetProperties()[1])

        # create & apply dirichlet bcs
        bcs_dirichlet_all = mp.CreateSubModelPart("BoundaryCondtionsDirichletAll")
        bcs_dirichlet_all.AddNodes([2,4])

        bcs_dirichlet_mv = mp.CreateSubModelPart("BoundaryCondtionsDirichletMove")
        bcs_dirichlet_mv.AddNodes([1,3])

        self._apply_dirichlet_BCs(bcs_dirichlet_all)
        self._apply_dirichlet_BCs(bcs_dirichlet_mv,fix_type='YZ')

        # create & apply neumann bcs
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[1],mp.GetProperties()[1])
        mp.CreateNewCondition("PointLoadCondition3D1N",2,[3],mp.GetProperties()[1])

        bcs_neumann = mp.CreateSubModelPart("BoundaryCondtionsNeumann")
        bcs_neumann.AddNodes([1,3])
        bcs_neumann.AddConditions([1,2])

        KratosMultiphysics.VariableUtils().SetScalarVar(StructuralMechanicsApplication.POINT_LOAD_X, 1000000.0, bcs_neumann.Nodes)

        # solve
        self._solve_static(mp)

        # check results
        self.assertAlmostEqual(mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), 0.58054148514004470,4)
        self.assertAlmostEqual(mp.Nodes[3].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), 0.15072065295319598,4)


class DynamicPatchTestMembrane(BasePatchTestMembrane):

    def test_membrane_3d3n_dynamic(self):

        displacement_results = [-0.007089438412033325, -0.02975285777563795,
         -0.06790344136536157, -0.11141464840612207, -0.14087911509274914,
          -0.13697941464148805, -0.10010593388219215, -0.051581204558573825,
           -0.013994666071832161, 0.0007937785999259335, -0.008609705698025766]


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

    def test_membrane_3d4n_dynamic(self):

        displacement_results = [-0.007965444599683935, -0.03429569272708872,
         -0.08347269985836575, -0.15492404513176913, -0.23652735237940312,
          -0.30707829314392593, -0.3272427614896494, -0.30090356124783313,
           -0.2666017458572726, -0.2435556601417446, -0.22710935966385354]


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


    def test_membrane_3d3n_dynamic_explicit(self):

        displacement_results = [-0.00040479957654093675, -0.0016324179587680046,
         -0.0037704865143040603, -0.006897010585329448, -0.01105166110793577,
          -0.016228737864184672, -0.02238079840616669, -0.029428730428419043,
           -0.037274683949580445, -0.045814472093039196, -0.054946932064362665,
            -0.06457910325238626, -0.07462741105917634, -0.08501594513456154,
             -0.09567325302012067, -0.10652892812279724, -0.11751086332538496,
             -0.1285435789039004, -0.13954765895506027, -0.1504401029901513,
              -0.16113531138791007, -0.17154643332610123, -0.18158686573505992,
               -0.19117176426047214, -0.20021949008993895, -0.20865296062028277,
                -0.21640089621275393, -0.22339896268355758, -0.22959080465276283,
                 -0.2349289543566711, -0.23937559005493414, -0.24290311265638176,
                  -0.24549451123781657, -0.24714349746297454, -0.2478544028992154,
                   -0.24764184823170488, -0.24653020619437496, -0.2445528888815245,
                    -0.24175149473282598, -0.23817485172406405, -0.23387799233376924,
                     -0.22892109368859137, -0.22336841348997402, -0.21728724908181823,
                      -0.21074694331627486, -0.20381795666297683, -0.19657102031530843,
                       -0.18907638002801846, -0.18140313534463032, -0.1736186740930271]


        current_model = KratosMultiphysics.Model()
        mp = self._set_up_system_3d3n(current_model,explicit_dynamics=True)


        #time integration parameters
        dt = 0.01
        time = 0.0
        end_time = 0.5
        step = 0

        self._set_and_fill_buffer(mp,3,dt)
        strategy_expl = _create_dynamic_explicit_strategy(mp,'central_differences')
        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)
            strategy_expl.Solve()
            self._check_dynamic_results(mp.Nodes[10],step-1,displacement_results)

def _create_dynamic_explicit_strategy(mp,scheme_name):
    if (scheme_name=='central_differences'):
        scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(0.00,0.00,0.00)
    elif scheme_name=='multi_stage':
        scheme = StructuralMechanicsApplication.ExplicitMultiStageKimScheme(0.33333333333333333)

    strategy = StructuralMechanicsApplication.MechanicalExplicitStrategy(mp,scheme,0,0,1)
    strategy.SetEchoLevel(0)
    return strategy

if __name__ == '__main__':
    KratosUnittest.main()
