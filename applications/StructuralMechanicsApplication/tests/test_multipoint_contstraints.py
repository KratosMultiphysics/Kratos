from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestMultipointConstraints(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _add_variables(self, mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

    def _add_dofs(self, mp):
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,
            mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,
            mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,
            mp)

        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.VELOCITY_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.VELOCITY_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.VELOCITY_Z, mp)

        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.ACCELERATION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.ACCELERATION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.ACCELERATION_Z, mp)

    def _apply_material_properties(self, mp, dim):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS, 210e9)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS, 1.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY, 1.0)

        g = [0, 0, 0]
        mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,
                                       g)

        if dim == 2:
            cl = KratosMultiphysics.StructuralMechanicsApplication.LinearElasticPlaneStrain2DLaw(
            )
        else:
            cl = KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3DLaw(
            )
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, cl)

    def _apply_BCs(self, mp):
        bcs = mp.GetSubModelPart("FixedEdgeNodes")
        KratosMultiphysics.VariableUtils().SetVariable(
            KratosMultiphysics.DISPLACEMENT_X, 0.0, bcs.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(
            KratosMultiphysics.DISPLACEMENT_Y, 0.0, bcs.Nodes)

        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_X, True, bcs.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_Y, True, bcs.Nodes)

        bcmn = mp.GetSubModelPart("MovingNodes")
        KratosMultiphysics.VariableUtils().SetVariable(
            KratosMultiphysics.DISPLACEMENT_X, 0.01, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(
            KratosMultiphysics.DISPLACEMENT_Y, 0.0, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_X, True, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_Y, True, bcmn.Nodes)

    def _setup_solver(self, mp):

        #define a minimal newton raphson solver
        self.linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        self.builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        self.scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(
            -0.01)
        self.convergence_criterion = KratosMultiphysics.ResidualCriteria(
            1e-10, 1e-12)
        self.convergence_criterion.SetEchoLevel(0)

        max_iters = 1000
        compute_reactions = False
        reform_step_dofs = True
        move_mesh_flag = False
        self.strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
            mp, self.scheme, self.linear_solver, self.convergence_criterion,
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

    def _check_results(self, mp):
        disp1 = mp.Nodes[16].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X, 0)
        disp2 = mp.Nodes[6].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(disp1, disp2, 5)

        disp1 = mp.Nodes[16].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y, 0)
        disp2 = mp.Nodes[6].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(disp1, disp2, 5)

        disp1 = 0.5 * (
            mp.Nodes[16].GetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT_X, 0) + mp.Nodes[18]
            .GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0))
        disp2 = mp.Nodes[7].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(disp1, disp2, 5)

        disp1 = 0.5 * (
            mp.Nodes[16].GetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT_Y, 0) + mp.Nodes[18]
            .GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0))
        disp2 = mp.Nodes[7].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(disp1, disp2, 5)

        disp1 = mp.Nodes[18].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X, 0)
        disp2 = mp.Nodes[9].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X, 0)
        self.assertAlmostEqual(disp1, disp2, 5)

        disp1 = mp.Nodes[18].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y, 0)
        disp2 = mp.Nodes[9].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y, 0)
        self.assertAlmostEqual(disp1, disp2, 5)

    def _setup_model_part(self, mp):
        #create nodes
        mp.CreateNewNode(1, 0.00000, 1.00000, 0.00000)
        mp.CreateNewNode(2, 0.00000, 0.50000, 0.00000)
        mp.CreateNewNode(3, 0.50000, 1.00000, 0.00000)
        mp.CreateNewNode(4, 0.50000, 0.50000, 0.00000)
        mp.CreateNewNode(5, 0.00000, 0.00000, 0.00000)
        mp.CreateNewNode(6, 1.00000, 1.00000, 0.00000)
        mp.CreateNewNode(7, 1.00000, 0.50000, 0.00000)
        mp.CreateNewNode(8, 0.50000, 0.00000, 0.00000)
        mp.CreateNewNode(9, 1.00000, 0.00000, 0.00000)
        mp.CreateNewNode(10, 1.50000, 1.00000, 0.00000)
        mp.CreateNewNode(12, 1.50000, 0.00000, 0.00000)
        mp.CreateNewNode(13, 2.00000, 1.00000, 0.00000)
        mp.CreateNewNode(15, 2.00000, 0.00000, 0.00000)
        mp.CreateNewNode(16, 1.00000, 1.00000, 0.00000)
        mp.CreateNewNode(18, 1.00000, 0.00000, 0.00000)

        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("FixedEdgeNodes")
        bcs.AddNodes([1, 2, 5, 13, 15])

        bcmn = mp.CreateSubModelPart("MovingNodes")
        bcmn.AddNodes([10, 12])

        #create Element
        mp.CreateNewElement("SmallDisplacementElement2D4N", 2,
                            [13, 10, 12, 15], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D4N", 4,
                            [10, 16, 18, 12], mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D4N", 5, [2, 4, 3, 1],
                            mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D4N", 6, [5, 8, 4, 2],
                            mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D4N", 7, [4, 7, 6, 3],
                            mp.GetProperties()[1])
        mp.CreateNewElement("SmallDisplacementElement2D4N", 8, [8, 9, 7, 4],
                            mp.GetProperties()[1])

    def _apply_mpc_constraints(self, mp, cm):
        cm.ExecuteAfterOutputStep()
        cm.AddMasterSlaveRelation(
            mp.Nodes[16], KratosMultiphysics.DISPLACEMENT_Y, mp.Nodes[6],
            KratosMultiphysics.DISPLACEMENT_Y, 1.0, 0)
        cm.AddMasterSlaveRelation(
            mp.Nodes[16], KratosMultiphysics.DISPLACEMENT_X, mp.Nodes[6],
            KratosMultiphysics.DISPLACEMENT_X, 1.0, 0)

        cm.AddMasterSlaveRelation(
            mp.Nodes[16], KratosMultiphysics.DISPLACEMENT_X, mp.Nodes[7],
            KratosMultiphysics.DISPLACEMENT_X, 0.5, 0)
        cm.AddMasterSlaveRelation(
            mp.Nodes[18], KratosMultiphysics.DISPLACEMENT_X, mp.Nodes[7],
            KratosMultiphysics.DISPLACEMENT_X, 0.5, 0)
        cm.AddMasterSlaveRelation(
            mp.Nodes[16], KratosMultiphysics.DISPLACEMENT_Y, mp.Nodes[7],
            KratosMultiphysics.DISPLACEMENT_Y, 0.5, 0)
        cm.AddMasterSlaveRelation(
            mp.Nodes[18], KratosMultiphysics.DISPLACEMENT_Y, mp.Nodes[7],
            KratosMultiphysics.DISPLACEMENT_Y, 0.5, 0)

        cm.AddMasterSlaveRelation(
            mp.Nodes[18], KratosMultiphysics.DISPLACEMENT_X, mp.Nodes[9],
            KratosMultiphysics.DISPLACEMENT_X, 1.0, 0)
        cm.AddMasterSlaveRelation(
            mp.Nodes[18], KratosMultiphysics.DISPLACEMENT_Y, mp.Nodes[9],
            KratosMultiphysics.DISPLACEMENT_Y, 1.0, 0)

    def _set_and_fill_buffer(self, mp, buffer_size, delta_time):
        # Set buffer size
        mp.SetBufferSize(buffer_size)
        # Fill buffer
        time = mp.ProcessInfo[KratosMultiphysics.TIME]
        time = time - delta_time * (buffer_size)
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for size in range(0, buffer_size):
            step = size - (buffer_size - 1)
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            mp.CloneTimeStep(time)

        mp.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    def test_MPC_Constraints(self):
        dim = 2
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")

        self._add_variables(mp)
        self._setup_model_part(mp)
        self._add_dofs(mp)
        self._apply_material_properties(mp, dim)

        #time integration parameters
        dt = 0.005
        time = 0.0
        end_time = 0.025
        step = 0

        self._set_and_fill_buffer(mp, 2, dt)
        # Applying boundary conditions
        self._apply_BCs(mp)
        # Applying constraints
        cm = KratosMultiphysics.StructuralMechanicsApplication.ApplyMultipointConstraintsProcess(
            mp)
        self._apply_mpc_constraints(mp, cm)
        # Solving the system of equations        
        self._setup_solver(mp)

        while (time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)
            if (step > 3):
                self._solve()
        # Checking the results
        self._check_results(mp)
        self._reset()


if __name__ == '__main__':
    KratosUnittest.main()
