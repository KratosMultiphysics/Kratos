from __future__ import print_function, absolute_import, division

import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, sin, cos, pi, exp, atan

class SpringDamperElementTests(KratosUnittest.TestCase):
    def setUp(self):
        pass
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)

    def _apply_material_properties(self,mp):
        cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _add_dofs(self,node):
        node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
        node.AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X)
        node.AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y)
        node.AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z)

    def _add_bcs(self,node):
        node.Fix(KratosMultiphysics.DISPLACEMENT_X)
        node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        node.Fix(KratosMultiphysics.ROTATION_X)
        node.Fix(KratosMultiphysics.ROTATION_Y)
        node.Fix(KratosMultiphysics.ROTATION_Z)

    def _apply_harmonic_cosine_load(self,amplitude,frequency,nodes,time):
        for node in nodes:
            force = amplitude * cos(2*pi*frequency*time)
            node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD,0,[0,force,0])

    def _apply_harmonic_sine_load(self,amplitude,frequency,nodes,time):
        for node in nodes:
            force = amplitude * sin(2*pi*frequency*time)
            node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD,0,[0,force,0])

    def _solve(self,mp):

        #define a minimal newton raphson dynamic solver
        damp_factor_m = -0.01
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m)
        # convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-14,1e-20)
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-4,1e-9)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 20
        compute_reactions = False
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

        strategy.Check()

        strategy.Solve()

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

    def _set_up_mdof_system(self, current_model):
        mp = current_model.CreateModelPart("mdof")
        self._add_variables(mp)
        self._apply_material_properties(mp)

        #create nodes
        n01 = mp.CreateNewNode(1,0.0,20.0,0.0)
        n02 = mp.CreateNewNode(2,0.0,10.0,0.0)
        n03 = mp.CreateNewNode(3,0.0,0.0,0.0)

        #create elements
        e01 = mp.CreateNewElement("NodalConcentratedElement3D1N",1, [1],  mp.GetProperties()[1])
        e02 = mp.CreateNewElement("SpringDamperElement3D2N",     2, [1,2],mp.GetProperties()[1])
        e03 = mp.CreateNewElement("NodalConcentratedElement3D1N",3, [2],  mp.GetProperties()[1])
        e04 = mp.CreateNewElement("SpringDamperElement3D2N",     4, [2,3],mp.GetProperties()[1])
        e05 = mp.CreateNewElement("SpringDamperElement3D2N",     5, [1,3],mp.GetProperties()[1])

        e01.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,0.0,0.0])
        e01.SetValue(KratosMultiphysics.NODAL_MASS,0.0)
        e02.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,0.0,0.0])
        e02.SetValue(KratosMultiphysics.NODAL_MASS,0.0)
        e03.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,0.0,0.0])
        e03.SetValue(KratosMultiphysics.NODAL_MASS,0.0)
        e04.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,0.0,0.0])
        e04.SetValue(KratosMultiphysics.NODAL_MASS,0.0)
        e05.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,0.0,0.0])
        e05.SetValue(KratosMultiphysics.NODAL_MASS,0.0)

        #add dofs and set bcs
        for node in mp.Nodes:
            self._add_dofs(node)
            self._add_bcs(node)
        n03.Fix(KratosMultiphysics.DISPLACEMENT_Y)

        return mp

    def _set_up_sdof_system(self, current_model):
        mp = current_model.CreateModelPart("sdof")
        self._add_variables(mp)
        self._apply_material_properties(mp)

        #create nodes
        n01 = mp.CreateNewNode(1,0.0,10.0,0.0)
        n02 = mp.CreateNewNode(2,0.0,0.0,0.0)

        #create elements
        e01 = mp.CreateNewElement("NodalConcentratedElement3D1N",1,[1],mp.GetProperties()[1])
        e02 = mp.CreateNewElement("SpringDamperElement3D2N",2,[1,2],mp.GetProperties()[1])

        e01.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,0.0,0.0])
        e01.SetValue(KratosMultiphysics.NODAL_MASS,0.0)
        e02.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,0.0,0.0])
        e02.SetValue(KratosMultiphysics.NODAL_MASS,0.0)

        #create condition
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[1],mp.GetProperties()[1])

        #add dofs and set bcs
        for node in mp.Nodes:
            self._add_dofs(node)
            self._add_bcs(node)
        n02.Fix(KratosMultiphysics.DISPLACEMENT_Y)

        return mp


    def test_undamped_mdof_system_dynamic(self):
        current_model = KratosMultiphysics.Model()
        mp = self._set_up_mdof_system(current_model)

        #set parameters
        mp.Elements[1].SetValue(KratosMultiphysics.NODAL_MASS,80.0)
        mp.Elements[2].SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,200.0,0.0])
        mp.Elements[3].SetValue(KratosMultiphysics.NODAL_MASS,8.0)
        mp.Elements[4].SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,125.0,0.0])

        #set initial conditions
        mp.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,-1.0)
        mp.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,1.0)

        #time integration parameters
        dt = 0.001
        time = 0.0
        end_time = 4.0
        step = 0

        self._set_and_fill_buffer(mp,2,dt)

        #parameters for analytical solution
        phi11 = 1.0
        phi12 = 1.0
        phi21 = 0.63
        phi22 = -15.88
        K11 = -0.901
        K12 = 0
        K21 = -0.099
        K22 = 0
        omega_E_1 = sqrt(0.93)
        omega_E_2 = sqrt(42.2)

        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)

            self._solve(mp)

            current_analytical_displacement_y_1 = phi11 * (K11*cos(omega_E_1*time) + K12*sin(omega_E_1*time)) \
                + phi12 * (K21*cos(omega_E_2*time) + K22*sin(omega_E_2*time))
            self.assertAlmostEqual(mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0),current_analytical_displacement_y_1,delta=1e-2)
            current_analytical_displacement_y_2 = phi21 * (K11*cos(omega_E_1*time) + K12*sin(omega_E_1*time)) \
                + phi22 * (K21*cos(omega_E_2*time) + K22*sin(omega_E_2*time))
            self.assertAlmostEqual(mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0),current_analytical_displacement_y_2,delta=1e-2)

    def test_undamped_sdof_system_harmonic(self):
        current_model = KratosMultiphysics.Model()
        mp = self._set_up_sdof_system(current_model)

        #set parameters
        mass = 80.0
        stiffness = 200.0
        F_0 = 100.0
        freq_ratio = 0.5 #beta
        eigen_freq = sqrt(stiffness/mass)  #omega
        excitation_freq = freq_ratio * eigen_freq #Omega
        init_u = 0.0
        init_v = 0.0
        mp.Elements[1].SetValue(KratosMultiphysics.NODAL_MASS,mass)
        mp.Elements[2].SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,stiffness,0.0])

        #time integration parameters
        dt = 0.01
        time = 0.0
        end_time = 10.0
        step = 0

        self._set_and_fill_buffer(mp,2,dt)

        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)
            self._apply_harmonic_sine_load(F_0,excitation_freq/(2*pi),[mp.Nodes[1]],time)

            self._solve(mp)

            u_trans = init_u * cos(eigen_freq*time) + init_v/eigen_freq * sin(eigen_freq*time) \
                - F_0/stiffness * freq_ratio / (1-freq_ratio*freq_ratio) * sin(eigen_freq*time)
            u_steadystate = F_0 / stiffness * 1 / (1-freq_ratio*freq_ratio) * sin(excitation_freq*time)
            current_analytical_displacement_y = u_trans + u_steadystate

            self.assertAlmostEqual(mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0), \
                current_analytical_displacement_y,delta=5e-3)


    def test_damped_mdof_system_dynamic(self):
        current_model = KratosMultiphysics.Model()

        mp = self._set_up_mdof_system(current_model)

        #set parameters
        mp.Elements[1].SetValue(KratosMultiphysics.NODAL_MASS,1.0)
        mp.Elements[2].SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,1.0,0.0])
        mp.Elements[2].SetValue(StructuralMechanicsApplication.NODAL_DAMPING_RATIO,[0.0,0.05,0.0])
        mp.Elements[3].SetValue(KratosMultiphysics.NODAL_MASS,2.0)
        mp.Elements[4].SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[0.0,2.0,0.0])
        mp.Elements[4].SetValue(StructuralMechanicsApplication.NODAL_DAMPING_RATIO,[0.0,0.4,0.0])
        mp.Elements[5].SetValue(StructuralMechanicsApplication.NODAL_DAMPING_RATIO,[0.0,0.15,0.0])

        #set initial conditions
        mp.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,2.0)
        mp.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,1.0)

        #time integration parameters
        dt = 0.05
        time = 0.0
        end_time = 20.0
        step = 0

        self._set_and_fill_buffer(mp,2,dt)

        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)

            self._solve(mp)
            current_analytical_displacement_y_1 = 2.0148 * exp(-0.08334*time) * sin(0.70221*time+atan(8.153))
            self.assertAlmostEqual(mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0), \
                current_analytical_displacement_y_1,delta=5e-2)

            current_analytical_displacement_y_2 = 1.0064 * exp(-0.08334*time) * sin(0.70221*time+atan(9.031))
            self.assertAlmostEqual(mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0), \
                current_analytical_displacement_y_2,delta=5e-2)

    def test_undamped_mdof_system_eigen(self):
        import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
        # FEAST is available otherwise this test is not being called
        if not hasattr(KratosMultiphysics.ExternalSolversApplication, "PastixSolver"):
            self.skipTest("Pastix Solver is not available")

        current_model = KratosMultiphysics.Model()

        mp = self._set_up_mdof_system(current_model)

        #set parameters
        mp.Elements[1].SetValue(KratosMultiphysics.NODAL_MASS,20.0)
        mp.Elements[2].SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[200.0,200.0,200.0])
        mp.Elements[3].SetValue(KratosMultiphysics.NODAL_MASS,40.0)
        mp.Elements[4].SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS,[400.0,400.0,400.0])

        #create solver
        eigen_solver_parameters = KratosMultiphysics.Parameters("""
            {
                "solver_type": "FEAST",
                "print_feast_output": false,
                "perform_stochastic_estimate": false,
                "solve_eigenvalue_problem": true,
                "lambda_min": 0.0,
                "lambda_max": 4.0e5,
                "number_of_eigenvalues": 2,
                "search_dimension": 18,
                "linear_solver_settings": {
                    "solver_type" : "pastix",
                    "echo_level" : 0
                }
            }""")
        feast_system_solver = ExternalSolversApplication.PastixComplexSolver(eigen_solver_parameters["linear_solver_settings"])
        eigen_solver = ExternalSolversApplication.FEASTSolver(eigen_solver_parameters, feast_system_solver)
        scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(eigen_solver)

        solver = StructuralMechanicsApplication.EigensolverStrategy(
            mp,
            scheme,
            builder_and_solver)

        solver.Solve()

        current_eigenvalues = [ev for ev in mp.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]]
        analytical_eigenvalues = [5,20]
        for ev in range(len(analytical_eigenvalues)):
            self.assertAlmostEqual(current_eigenvalues[ev], analytical_eigenvalues[ev])

if __name__ == '__main__':
    KratosUnittest.main()
