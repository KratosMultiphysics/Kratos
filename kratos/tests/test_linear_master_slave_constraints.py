from __future__ import print_function, absolute_import, division

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils

structural_mechanics_is_available = KratosUtils.CheckIfApplicationsAvailable("StructuralMechanicsApplication")
if structural_mechanics_is_available:
    import KratosMultiphysics.StructuralMechanicsApplication as SMA


class TestLinearMultipointConstraints(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_variables(self):
        self.mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.mp.AddNodalSolutionStepVariable(KM.REACTION)
        self.mp.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.mp.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.mp.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)

    def _add_dofs(self):
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, self.mp)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, self.mp)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, self.mp)

        KM.VariableUtils().AddDof(KM.VELOCITY_X, self.mp)
        KM.VariableUtils().AddDof(KM.VELOCITY_Y, self.mp)
        KM.VariableUtils().AddDof(KM.VELOCITY_Z, self.mp)

        KM.VariableUtils().AddDof(KM.ACCELERATION_X, self.mp)
        KM.VariableUtils().AddDof(KM.ACCELERATION_Y, self.mp)
        KM.VariableUtils().AddDof(KM.ACCELERATION_Z, self.mp)

    def _apply_material_properties(self, dim):
        #define properties
        self.mp.GetProperties()[1].SetValue(KM.YOUNG_MODULUS, 210e9)
        self.mp.GetProperties()[1].SetValue(KM.POISSON_RATIO, 0.3)
        self.mp.GetProperties()[1].SetValue(KM.THICKNESS, 1.0)
        self.mp.GetProperties()[1].SetValue(KM.DENSITY, 1.0)

        g = [0, 0, 0]
        self.mp.GetProperties()[1].SetValue(KM.VOLUME_ACCELERATION, g)

        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = dim
        if dim == 2:
            cl = KM.StructuralMechanicsApplication.LinearElasticPlaneStrain2DLaw()
        else:
            cl = KM.StructuralMechanicsApplication.LinearElastic3DLaw()
        self.mp.GetProperties()[1].SetValue(KM.CONSTITUTIVE_LAW, cl)

    def _apply_BCs(self):
        bcs = self.mp.GetSubModelPart("FixedEdgeNodes")
        KM.VariableUtils().SetScalarVar(KM.DISPLACEMENT_X, 0.0, bcs.Nodes)
        KM.VariableUtils().SetScalarVar(KM.DISPLACEMENT_Y, 0.0, bcs.Nodes)

        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_X, True, bcs.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Y, True, bcs.Nodes)

        bcmn = self.mp.GetSubModelPart("MovingNodes")
        KM.VariableUtils().SetScalarVar(KM.DISPLACEMENT_X, 0.01, bcmn.Nodes)
        KM.VariableUtils().SetScalarVar(KM.DISPLACEMENT_Y, 0.0, bcmn.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_X, True, bcmn.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Y, True, bcmn.Nodes)

    def _setup_solver(self, solving_with = "Block", linear_solver = "AMGCL"):

        #define a minimal newton raphson solver
        if linear_solver == "AMGCL":
            self.linear_solver = KM.AMGCLSolver()
        else:
            self.linear_solver = KM.SkylineLUFactorizationSolver()
        if solving_with == "Block":
            self.builder_and_solver = KM.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        elif solving_with == "LM":
            self.builder_and_solver = KM.ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier(self.linear_solver, True, False, False)
        elif solving_with == "DoubleLM":
            self.builder_and_solver = KM.ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier(self.linear_solver, True, False, True)
        else: # Block default
            self.builder_and_solver = KM.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        self.scheme = KM.ResidualBasedBossakDisplacementScheme(-0.01)
        self.convergence_criterion = KM.ResidualCriteria(1e-6, 1e-9)
        self.convergence_criterion.SetEchoLevel(0)

        max_iters = 10
        compute_reactions = True
        reform_step_dofs = True
        move_mesh_flag = False
        self.strategy = KM.ResidualBasedNewtonRaphsonStrategy(self.mp, self.scheme, self.linear_solver, self.convergence_criterion, self.builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag)
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

    def _basic_check_results(self):
        reactionx1 = self.mp.Nodes[1].GetSolutionStepValue(KM.REACTION_X, 0)
        self.assertLessEqual(abs((reactionx1 - -1413464323.8223937)/(-1413464323.8223937)), 1.0e-2)
        reactiony1 = self.mp.Nodes[1].GetSolutionStepValue(KM.REACTION_Y, 0)
        self.assertLessEqual(abs((reactiony1 - -605769230.7692306)/(-605769230.7692306)), 1.0e-2)

        reactionx4 = self.mp.Nodes[4].GetSolutionStepValue(KM.REACTION_X, 0)
        self.assertLessEqual(abs((reactionx4 - -1413467109.1832492)/(-1413467109.1832492)), 1.0e-2)
        reactiony4 = self.mp.Nodes[4].GetSolutionStepValue(KM.REACTION_Y, 0)
        self.assertLessEqual(abs((reactiony4 - 605769230.7692306)/(605769230.7692306)), 1.0e-2)

        dispx3 = self.mp.Nodes[3].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx3, 0.01, 4)
        dispy3 = self.mp.Nodes[3].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy3, 0.0, 4)

        dispx2 = self.mp.Nodes[2].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx2, 0.01, 4)
        dispy2 = self.mp.Nodes[2].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy2, 0.0, 4)

        dispx3 = self.mp.Nodes[3].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx3, 0.01, 4)
        dispy3 = self.mp.Nodes[3].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy3, 0.0, 4)

    def _advanced_check_results(self):
        dispx13 = self.mp.Nodes[13].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx13, 0.01, 4)
        dispy13 = self.mp.Nodes[13].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy13, 0.0, 4)

        dispx14 = self.mp.Nodes[14].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx14, 0.01, 4)
        dispy14 = self.mp.Nodes[14].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy14, 0.0, 4)

        dispx15 = self.mp.Nodes[15].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx15, 0.01, 4)
        dispy15 = self.mp.Nodes[15].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy15, 0.0, 4)

        dispx11 = self.mp.Nodes[11].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx11, 0.0077238, 4)
        dispy11 = self.mp.Nodes[11].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy11, 0.0, 4)

        dispx4 = self.mp.Nodes[4].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(dispx4, 0.0022754, 4)
        dispy4 = self.mp.Nodes[4].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(dispy4, 0.0, 4)


        disp1 = self.mp.Nodes[16].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        disp2 = self.mp.Nodes[6].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(disp1, disp2, 4)
        self.assertAlmostEqual(disp1, 0.0049994, 4)
        self.assertAlmostEqual(disp2, 0.0049994, 4)
        #print("Test 1 :: ", disp1," == ",disp2)

        disp1 = self.mp.Nodes[16].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        disp2 = self.mp.Nodes[6].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(disp1, disp2, 5)
        self.assertAlmostEqual(disp1, -0.0011584, 4)
        self.assertAlmostEqual(disp1, -0.0011584, 4)
        #print("Test 2 :: ", disp1," == ",disp2)

        disp1 = self.mp.Nodes[7].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        disp2 = self.mp.Nodes[17].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(disp1, disp2, 4)
        self.assertAlmostEqual(disp1, 0.0049994, 4)
        self.assertAlmostEqual(disp2, 0.0049994, 4)
        #print("Test 3 :: ", disp1," == ",disp2)

        disp1 = self.mp.Nodes[7].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        disp2 = self.mp.Nodes[17].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(disp1, disp2, 4)
        #print("Test 4 :: ", disp1," == ",disp2)

        disp1 = self.mp.Nodes[18].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        disp2 = self.mp.Nodes[9].GetSolutionStepValue(KM.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(disp1, disp2, 4)
        self.assertAlmostEqual(disp1, 0.0049994, 4)
        self.assertAlmostEqual(disp1, 0.0049994, 4)
        #print("Test 5 :: ", disp1," == ",disp2)

        disp1 = self.mp.Nodes[18].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        disp2 = self.mp.Nodes[9].GetSolutionStepValue(KM.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(disp1, disp2, 4)
        self.assertAlmostEqual(disp1, 0.0011584, 4)
        self.assertAlmostEqual(disp1, 0.0011584, 4)
        #print("Test 6 :: ", disp1," == ",disp2)

    def _basic_setup_model_part(self):
        #create nodes
        self.mp.CreateNewNode(1, 0.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(2, 1.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(3, 1.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(4, 0.00000, 1.00000, 0.00000)

        #create a submodelpart for boundary conditions
        bcs = self.mp.CreateSubModelPart("FixedEdgeNodes")
        bcs.AddNodes([1, 4])

        bcmn = self.mp.CreateSubModelPart("MovingNodes")
        bcmn.AddNodes([2])

        #create Element
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 1, [1, 2, 3, 4], self.mp.GetProperties()[1])

    def _advanced_setup_model_part(self):
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
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 1, [14, 11, 12, 15], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 2, [13, 10, 11, 14], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 3, [11, 17, 18, 12], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 4, [10, 16, 17, 11], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 5, [2, 4, 3, 1], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 6, [5, 8, 4, 2], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 7, [4, 7, 6, 3], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("SmallDisplacementElement2D4N", 8, [8, 9, 7, 4], self.mp.GetProperties()[1])

    def _basic_apply_mpc_constraints(self):
        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, self.mp.Nodes[2], KM.DISPLACEMENT_X, self.mp.Nodes[3], KM.DISPLACEMENT_X, 1.0, 0)
        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, self.mp.Nodes[2], KM.DISPLACEMENT_Y, self.mp.Nodes[3], KM.DISPLACEMENT_Y, 1.0, 0)

    def _advanced_apply_mpc_constraints(self):

        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, self.mp.Nodes[16], KM.DISPLACEMENT_X, self.mp.Nodes[6], KM.DISPLACEMENT_X, 1.0, 0)
        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, self.mp.Nodes[16], KM.DISPLACEMENT_Y, self.mp.Nodes[6], KM.DISPLACEMENT_Y, 1.0, 0)


        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3, self.mp.Nodes[17], KM.DISPLACEMENT_X, self.mp.Nodes[7], KM.DISPLACEMENT_X, 1.0, 0)
        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 4, self.mp.Nodes[17], KM.DISPLACEMENT_Y, self.mp.Nodes[7], KM.DISPLACEMENT_Y, 1.0, 0)

        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 5, self.mp.Nodes[18], KM.DISPLACEMENT_X, self.mp.Nodes[9], KM.DISPLACEMENT_X, 1.0, 0)
        self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 6, self.mp.Nodes[18], KM.DISPLACEMENT_Y, self.mp.Nodes[9], KM.DISPLACEMENT_Y, 1.0, 0)

    def _set_and_fill_buffer(self, buffer_size, delta_time):
        # Set buffer size
        self.mp.SetBufferSize(buffer_size)
        # Fill buffer
        time = self.mp.ProcessInfo[KM.TIME]
        time = time - delta_time * (buffer_size)
        self.mp.ProcessInfo.SetValue(KM.TIME, time)
        for size in range(0, buffer_size):
            step = size - (buffer_size - 1)
            self.mp.ProcessInfo.SetValue(KM.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            self.mp.CloneTimeStep(time)

        self.mp.ProcessInfo[KM.IS_RESTARTED] = False

    def _basic_setup_test(self, solving_with = "Block", linear_solver = "AMGCL"):
        dim = 2
        current_model = KM.Model()
        self.mp= current_model.CreateModelPart("MainModelPart")
        self._add_variables()
        self._basic_setup_model_part()
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
        self._basic_apply_mpc_constraints()
        # Solving the system of equations
        self._setup_solver(solving_with, linear_solver)

        while (time <= end_time):
            time = time + dt
            step = step + 1
            self.mp.CloneTimeStep(time)
            self._solve()
        # Checking the results
        self._basic_check_results()
        self._reset()

    def _advanced_setup_test(self, solving_with = "Block", linear_solver = "AMGCL"):
        dim = 2
        current_model = KM.Model()
        self.mp= current_model.CreateModelPart("MainModelPart")
        self._add_variables()
        self._advanced_setup_model_part()
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
        self._advanced_apply_mpc_constraints()
        # Solving the system of equations
        self._setup_solver(solving_with, linear_solver)

        while (time <= end_time):
            time = time + dt
            step = step + 1
            self.mp.CloneTimeStep(time)
            self._solve()
        # Checking the results
        self._advanced_check_results()
        self._reset()

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_basic_MPC_Constraints(self):
        self._basic_setup_test("Block")

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_advanced_MPC_Constraints(self):
        self._advanced_setup_test("Block")

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_basic_LM_MPC_Constraints(self):
        self._basic_setup_test("LM", "LU")

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_advanced_LM_MPC_Constraints(self):
        self._advanced_setup_test("LM", "LU")

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_basic_Double_LM_MPC_Constraints(self):
        self._basic_setup_test("DoubleLM")

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_advanced_Double_LM_MPC_Constraints(self):
        self._advanced_setup_test("DoubleLM")

if __name__ == '__main__':
    KratosUnittest.main()
