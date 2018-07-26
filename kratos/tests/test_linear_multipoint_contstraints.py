from __future__ import print_function, absolute_import, division

# Importing the Kratos Library
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
try:
    import KratosMultiphysics.StructuralMechanicsApplication
    missing_external_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''',
                                    '{0}'.format(e)).group(1)

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
        KratosMultiphysics.VariableUtils().SetScalarVar(
            KratosMultiphysics.DISPLACEMENT_X, 0.0, bcs.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(
            KratosMultiphysics.DISPLACEMENT_Y, 0.0, bcs.Nodes)

        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_X, True, bcs.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_Y, True, bcs.Nodes)

        bcmn = self.mp.GetSubModelPart("MovingNodes")
        KratosMultiphysics.VariableUtils().SetScalarVar(
            KratosMultiphysics.DISPLACEMENT_X, 0.01, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(
            KratosMultiphysics.DISPLACEMENT_Y, 0.0, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_X, True, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_Y, True, bcmn.Nodes)

    def _setup_solver(self):

        #define a minimal newton raphson solver
        self.linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        #self.builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        self.builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolverWithConstraints(
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

    def test_MPC_Constraints(self):
        dim = 2
        self.mp = KratosMultiphysics.ModelPart("MainModelPart")
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


if __name__ == '__main__':
    KratosUnittest.main()
