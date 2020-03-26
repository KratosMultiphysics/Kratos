from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils
import math

structural_mechanics_is_available = KratosUtils.CheckIfApplicationsAvailable("StructuralMechanicsApplication")
if structural_mechanics_is_available:
    import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def inner_prod(a,b):
    tmp = 0
    for i in range(len(a)):
        tmp += a[i]*b[i]
    return tmp

class TestLinearMultipointConstraints(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_variables(self):
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

    def _add_dofs(self):
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,
            self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,
            self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,
            self.mp)

        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.VELOCITY_X, self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.VELOCITY_Y, self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.VELOCITY_Z, self.mp)

        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.ACCELERATION_X, self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.ACCELERATION_Y, self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.ACCELERATION_Z, self.mp)

    def _apply_material_properties(self, dim):
        #define properties
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS, 210e9)
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS, 1.0)
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY, 1.0)

        g = [0, 0, 0]
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,
                                       g)

        if dim == 2:
            cl = KratosMultiphysics.StructuralMechanicsApplication.LinearElasticPlaneStrain2DLaw(
            )
        else:
            cl = KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3DLaw(
            )
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, cl)

    def _apply_BCs(self):
        bcs = self.mp.GetSubModelPart("FixedEdgeNodes")
        KratosMultiphysics.VariableUtils().SetVariable(
            KratosMultiphysics.DISPLACEMENT_X, 0.0, bcs.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(
            KratosMultiphysics.DISPLACEMENT_Y, 0.0, bcs.Nodes)

        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_X, True, bcs.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_Y, True, bcs.Nodes)

        bcmn = self.mp.GetSubModelPart("MovingNodes")
        KratosMultiphysics.VariableUtils().SetVariable(
            KratosMultiphysics.DISPLACEMENT_X, 0.01, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(
            KratosMultiphysics.DISPLACEMENT_Y, 0.0, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_X, True, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_Y, True, bcmn.Nodes)

    def _setup_solver(self):

        #define a minimal newton raphson solver
        self.linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        self.builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(
            self.linear_solver)
        self.scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(
            -0.01)
        self.convergence_criterion = KratosMultiphysics.ResidualCriteria(
            1e-10, 1e-12)
        self.convergence_criterion.SetEchoLevel(0)

        max_iters = 100
        compute_reactions = False
        reform_step_dofs = True
        move_mesh_flag = False
        self.strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
            self.mp, self.scheme, self.linear_solver, self.convergence_criterion,
            self.builder_and_solver, max_iters, compute_reactions,
            reform_step_dofs, move_mesh_flag)
        self.strategy.SetEchoLevel(0)
        self.strategy.Initialize()

        self.strategy.Check()

    def _reset(self):
        del self.strategy
        del self.linear_solver
        del self.builder_and_solver
        del self.scheme
        del self.convergence_criterion

    def _solve(self):
        self.strategy.Solve()

    def _check_results(self):
        dispx13 = self.mp.Nodes[13].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx13, 0.01, 4)
        dispy13 = self.mp.Nodes[13].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy13, 0.0, 4)

        dispx14 = self.mp.Nodes[14].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx14, 0.01, 4)
        dispy14 = self.mp.Nodes[14].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy14, 0.0, 4)

        dispx15 = self.mp.Nodes[15].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx15, 0.01, 4)
        dispy15 = self.mp.Nodes[15].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy15, 0.0, 4)

        dispx11 = self.mp.Nodes[11].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx11, 0.0077238, 4)
        dispy11 = self.mp.Nodes[11].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy11, 0.0, 4)

        dispx4 = self.mp.Nodes[4].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx4, 0.0022754, 4)
        dispy4 = self.mp.Nodes[4].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy4, 0.0, 4)


        disp1 = self.mp.Nodes[16].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        disp2 = self.mp.Nodes[6].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(disp1, disp2, 4)
        self.assertAlmostEqual(disp1, 0.0049994, 4)
        self.assertAlmostEqual(disp2, 0.0049994, 4)
        #print("Test 1 :: ", disp1," == ",disp2)

        disp1 = self.mp.Nodes[16].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        disp2 = self.mp.Nodes[6].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(disp1, disp2, 5)
        self.assertAlmostEqual(disp1, -0.0011584, 4)
        self.assertAlmostEqual(disp1, -0.0011584, 4)
        #print("Test 2 :: ", disp1," == ",disp2)

        disp1 = self.mp.Nodes[7].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        disp2 = self.mp.Nodes[17].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(disp1, disp2, 4)
        self.assertAlmostEqual(disp1, 0.0049994, 4)
        self.assertAlmostEqual(disp2, 0.0049994, 4)
        #print("Test 3 :: ", disp1," == ",disp2)

        disp1 = self.mp.Nodes[7].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        disp2 = self.mp.Nodes[17].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(disp1, disp2, 4)
        #print("Test 4 :: ", disp1," == ",disp2)

        disp1 = self.mp.Nodes[18].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        disp2 = self.mp.Nodes[9].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(disp1, disp2, 4)
        self.assertAlmostEqual(disp1, 0.0049994, 4)
        self.assertAlmostEqual(disp1, 0.0049994, 4)
        #print("Test 5 :: ", disp1," == ",disp2)

        disp1 = self.mp.Nodes[18].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        disp2 = self.mp.Nodes[9].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(disp1, disp2, 4)
        self.assertAlmostEqual(disp1, 0.0011584, 4)
        self.assertAlmostEqual(disp1, 0.0011584, 4)
        #print("Test 6 :: ", disp1," == ",disp2)

    def _setup_model_part(self):
        #create nodes
        self.mp.CreateNewNode(1, 0.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(2, 0.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(3, 0.50000, 1.00000, 0.00000)
        self.mp.CreateNewNode(4, 0.50000, 0.50000, 0.00000)
        self.mp.CreateNewNode(5, 0.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(6, 1.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(7, 1.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(8, 0.50000, 0.00000, 0.00000)
        self.mp.CreateNewNode(9, 1.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(10, 1.50000, 1.00000, 0.00000)
        self.mp.CreateNewNode(11, 1.50000, 0.50000, 0.00000)
        self.mp.CreateNewNode(12, 1.50000, 0.00000, 0.00000)
        self.mp.CreateNewNode(13, 2.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(14, 2.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(15, 2.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(16, 1.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(17, 1.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(18, 1.00000, 0.00000, 0.00000)

        #create a submodelpart for boundary conditions
        bcs = self.mp.CreateSubModelPart("FixedEdgeNodes")
        bcs.AddNodes([1, 2, 5])

        bcmn = self.mp.CreateSubModelPart("MovingNodes")
        bcmn.AddNodes([13, 14, 15])

        #create Element
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 1,
                            [14, 11, 12, 15], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 2,
                            [13, 10, 11, 14], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 3,
                            [11, 17, 18, 12], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 4,
                            [10, 16, 17, 11], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 5,
                            [2, 4, 3, 1], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 6, [5, 8, 4, 2],
                            self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 7, [4, 7, 6, 3],
                            self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 8, [8, 9, 7, 4],
                            self.mp.GetProperties()[1])

    def _apply_mpc_constraints(self):

        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, self.mp.Nodes[16], KratosMultiphysics.DISPLACEMENT_X, self.mp.Nodes[6], KratosMultiphysics.DISPLACEMENT_X, 1.0, 0)
        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, self.mp.Nodes[16], KratosMultiphysics.DISPLACEMENT_Y, self.mp.Nodes[6], KratosMultiphysics.DISPLACEMENT_Y, 1.0, 0)


        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3, self.mp.Nodes[17], KratosMultiphysics.DISPLACEMENT_X, self.mp.Nodes[7], KratosMultiphysics.DISPLACEMENT_X, 1.0, 0)
        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 4, self.mp.Nodes[17], KratosMultiphysics.DISPLACEMENT_Y, self.mp.Nodes[7], KratosMultiphysics.DISPLACEMENT_Y, 1.0, 0)

        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 5, self.mp.Nodes[18], KratosMultiphysics.DISPLACEMENT_X, self.mp.Nodes[9], KratosMultiphysics.DISPLACEMENT_X, 1.0, 0)
        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 6, self.mp.Nodes[18], KratosMultiphysics.DISPLACEMENT_Y, self.mp.Nodes[9], KratosMultiphysics.DISPLACEMENT_Y, 1.0, 0)

    def _set_and_fill_buffer(self, buffer_size, delta_time):
        # Set buffer size
        self.mp.SetBufferSize(buffer_size)
        # Fill buffer
        time = self.mp.ProcessInfo[KratosMultiphysics.TIME]
        time = time - delta_time * (buffer_size)
        self.mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for size in range(0, buffer_size):
            step = size - (buffer_size - 1)
            self.mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            self.mp.CloneTimeStep(time)

        self.mp.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_MPC_Constraints(self):
        dim = 2
        current_model = KratosMultiphysics.Model()
        self.mp= current_model.CreateModelPart("MainModelPart")
        self.mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        self._add_variables()
        self._setup_model_part()
        self._add_dofs()
        self._apply_material_properties(dim)

        #time integration parameters
        dt = 0.002
        time = 0.0
        end_time = 0.01
        step = 0

        self._set_and_fill_buffer(2, dt)
        # Applying boundary conditions
        self._apply_BCs()
        # Applying constraints
        self._apply_mpc_constraints()
        # Solving the system of equations
        self._setup_solver()

        while (time <= end_time):
            time = time + dt
            step = step + 1
            self.mp.CloneTimeStep(time)
            self._solve()
        # Checking the results
        self._check_results()
        self._reset()


class TestLinearConstraints(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

    def _apply_material_properties(self, mp, dim, small_strain = True):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,210)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,1.0)

        g = [0,-10.0,0]
        mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = StructuralMechanicsApplication.LinearElasticPlaneStrain2DLaw()
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _set_buffer(self,mp):
        buffer_size = 3
        mp.SetBufferSize(buffer_size)
        # Cycle the buffer. This sets all historical nodal solution step data to
        # the current value and initializes the time stepping in the process info.
        mp.ProcessInfo[KratosMultiphysics.DELTA_TIME] = 1.0
        delta_time = mp.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        time = mp.ProcessInfo[KratosMultiphysics.TIME]
        step =-buffer_size
        time = time - delta_time * buffer_size
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for i in range(0, buffer_size):
            step = step + 1
            time = time + delta_time
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            mp.CloneTimeStep(time)

    def _apply_BCs(self,mp,A,b):
        pass

    def _create_strategy(self, mp):
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.DisplacementCriteria(1e-13,1e-14)
        convergence_criterion.SetEchoLevel(0)

        #max_iters = 1
        max_iters = 20
        compute_reactions = True
        reform_step_dofs = True
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

        return strategy

    def _solve_with_strategy(self, strategy, step):
        strategy.Check()
        strategy.Initialize()
        strategy.InitializeSolutionStep()
        strategy.Predict()
        strategy.SolveSolutionStep()
        strategy.FinalizeSolutionStep()

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_constraints(self):
        dim = 2

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("tl_solid_part")

        self._add_variables(mp)
        self._apply_material_properties(mp, dim, False)

        #
        #     3
        #  4     2
        #     1
        # Create nodes
        n1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        n2 = mp.CreateNewNode(2, 2.0, 1.0, 0.0)
        n3 = mp.CreateNewNode(3, 0.0, 2.0, 0.0)
        n4 = mp.CreateNewNode(4, -2.0, 1.0, 0.0)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        mp.CreateNewElement("TotalLagrangianElement2D4N", 1, [1,2,3,4], mp.GetProperties()[1])

        #node 1 is completely fixed
        n1.Fix(KratosMultiphysics.DISPLACEMENT_X)
        n1.Fix(KratosMultiphysics.DISPLACEMENT_Y)

        ################# apply a constraint on node 2
        #node2 is at allowed to slide normally to normal2 (taken as normal to the side 1-2)
        dx = n2.X - n1.X
        dy = n2.Y - n1.Y
        normal_2 = KratosMultiphysics.Array3([dy,-dx,0])/math.sqrt(dx**2+dy**2)
        inner_prod(normal_2,normal_2)
        master_dofs_2 = [n2.GetDof(KratosMultiphysics.DISPLACEMENT_Y)]
        slave_dofs_2 = [n2.GetDof(KratosMultiphysics.DISPLACEMENT_X)]
        RelationMatrix2 = KratosMultiphysics.Matrix(1,1)
        RelationMatrix2[0,0] = -normal_2[1]/normal_2[0]
        ConstantVector = KratosMultiphysics.Vector([0.0])

        constraint_2 = KratosMultiphysics.LinearMasterSlaveConstraint(2,

            master_dofs_2,
            slave_dofs_2,
            RelationMatrix2,
            ConstantVector)
        mp.AddMasterSlaveConstraint(constraint_2)

        ################# apply a slip constraint on node 4
        # note that SlipConstraint is practically a placeholder for the operations done just above
        #node4 is at allowed to slide normally to normal4
        dx = n1.X - n4.X
        dy = n1.Y - n4.Y
        normal_4 = KratosMultiphysics.Array3([dy,-dx,0])/math.sqrt(dx**2+dy**2)
        n4.SetValue(KratosMultiphysics.NORMAL, normal_4) #artificially set the normal

        constraint_4 = KratosMultiphysics.SlipConstraint(4,
            n4.GetDof(KratosMultiphysics.DISPLACEMENT_X),
            n4.GetDof(KratosMultiphysics.DISPLACEMENT_Y),
            n4.GetValue(KratosMultiphysics.NORMAL)
            )
        mp.AddMasterSlaveConstraint(constraint_4)

        #solve the problem
        strategy = self._create_strategy(mp)
        self._solve_with_strategy(strategy,0)

        ##verify the results
        d2 = n2.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        d3 = n3.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        d4 = n4.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)

        self.assertAlmostEqual(d2[0]*normal_2[0] + d2[1]*normal_2[1],0.0, 15)
        self.assertAlmostEqual(d4[0]*normal_4[0] + d4[1]*normal_4[1],0.0, 15)
        self.assertAlmostEqual(d3[0], 0.0, 15) #symmetry condition
        self.assertAlmostEqual(d3[1], 2.0*d2[1], 15)

        R2 = n2.GetSolutionStepValue(KratosMultiphysics.REACTION)
        R4 = n4.GetSolutionStepValue(KratosMultiphysics.REACTION)

        nR2 = inner_prod(normal_2, R2)
        nR4 = inner_prod(normal_4, R4)

        #check that Reactions are in compression
        self.assertTrue(nR2 < 0.0)
        self.assertTrue(nR4 < 0.0)

        #check that tangential component is zero
        tang_2 = R2 - nR2*normal_2
        tang_4 = R4 - nR4*normal_4

        self.assertTrue(tang_2.norm_2() < 1e-12)
        self.assertTrue(tang_4.norm_2() < 1e-12)

        self.assertEqual(mp.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER], 5) #4 if using Residual Criteria

if __name__ == '__main__':
    KratosUnittest.main()
