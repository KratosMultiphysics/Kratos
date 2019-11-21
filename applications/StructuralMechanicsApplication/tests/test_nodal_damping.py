from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, atan, cos, exp

class NodalDampingTests(KratosUnittest.TestCase):
    def setUp(self):
        pass
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)


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
        mp.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        for size in range(0, buffer_size):
            step = size - (buffer_size -1)
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            mp.CloneTimeStep(time)

        mp.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    def test_nodal_damping(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("sdof")

        self._add_variables(mp)

        #create node
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
        element = mp.CreateNewElement("NodalConcentratedDampedElement3D1N", 1, [1], None)
        mass = 1.0
        stiffness = 10.0
        damping = 1.0
        element.SetValue(KratosMultiphysics.NODAL_MASS,mass)
        element.SetValue(StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS,[0,stiffness,0])
        element.SetValue(StructuralMechanicsApplication.NODAL_DAMPING_RATIO,[0,damping,0])

        #time integration parameters
        dt = 0.005
        time = 0.0
        end_time = 5.0
        step = 0

        self._set_and_fill_buffer(mp,2,dt)

        #parameters for analytical solution
        omega = sqrt(stiffness/mass)
        D = damping / (2*mass*omega)
        omega_D = omega * sqrt(1-D*D)
        delta = damping / (2*mass)
        theta = atan(-(init_velocity+init_displacement*delta) / (omega_D*init_displacement))
        A = sqrt(init_displacement*init_displacement + ((init_velocity+init_displacement*delta) / omega_D)**2)

        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)

            self._solve(mp)
            current_analytical_displacement_y = A * cos(omega_D*time+theta) * exp(-delta*time)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0),current_analytical_displacement_y,delta=1e-3)


if __name__ == '__main__':
    KratosUnittest.main()
