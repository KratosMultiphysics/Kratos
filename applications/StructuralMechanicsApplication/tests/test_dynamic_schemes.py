from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, cos

class DynamicSchemesTests(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)        

        
    def _create_solver(self,mp, scheme_name):
        #define a minimal newton raphson dynamic solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        if (scheme_name == "newmark"):
            scheme = KratosMultiphysics.ResidualBasedNewmarkDisplacementScheme()
        elif (scheme_name == "bdf2"):
            scheme = KratosMultiphysics.ResidualBasedBDF2DisplacementScheme()
        else:
            damp_factor_m = -0.01
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
        
        return strategy

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
   
    def _base_test_dynamic_schemes(self, scheme_name = "bossak"):
        mp = KratosMultiphysics.ModelPart("sdof")
        self._add_variables(mp)

        # Create node
        node = mp.CreateNewNode(1,0.0,0.0,0.0)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)

        #add bcs and initial values
        init_displacement = 0.1
        init_velocity = 0.0
        node.Fix(KratosMultiphysics.DISPLACEMENT_X)
        node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,init_displacement)

        #create element
        element = mp.CreateNewElement("NodalConcentratedElement3D1N", 1, [1], mp.GetProperties()[1])
        mass = 1.0
        stiffness = 10.0
        damping = 1.0
        element.SetValue(KratosMultiphysics.NODAL_MASS, mass)
        element.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS, [0, stiffness,0])

        #time integration parameters
        dt = 0.01
        time = 0.0
        end_time = 0.1
        step = 0

        self._set_and_fill_buffer(mp,2,dt)

        #parameters for analytical solution
        omega = sqrt(stiffness/mass)
        A = init_displacement

        self.strategy = self._create_solver(mp, scheme_name)

        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)

            self.strategy.Solve()
            current_analytical_displacement_y = A * cos(omega*time)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0), current_analytical_displacement_y, delta=1e-3)
            
    def test_bossak_scheme(self):
        self._base_test_dynamic_schemes("bossak")
        
    def test_newmark_scheme(self):
        self._base_test_dynamic_schemes("newmark")
        
    #def test_bdf2_scheme(self):
        #self._base_test_dynamic_schemes("bdf2")

if __name__ == '__main__':
    KratosUnittest.main()
