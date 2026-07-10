import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, atan, cos, exp

class NodalDampingTests(KratosUnittest.TestCase):
    def setUp(self):
        pass
    def _add_variables(self,mp,explicit_dynamics=False):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
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
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        mp.CreateNewProperties(1)
        
        self._add_variables(mp)

        # Create node
        node = mp.CreateNewNode(1,0.0,0.0,0.0)
        _add_dofs(node)

        # Add bcs and initial values
        init_displacement = 0.1
        init_velocity = 0.0
        _set_dirichlet_bc(node,init_displacement,init_velocity)

        # Create element
        mass,stiffness,damping = _set_material_properties()
        _create_element(mp,mass,stiffness,damping)

        # Time integration parameters
        dt = 0.005
        time = 0.0
        end_time = 5.0
        step = 0

        self._set_and_fill_buffer(mp,2,dt)

        # Parameters for analytical solution
        omega_D,delta,theta,A = _return_parameters_analytical_solution(stiffness,mass,damping,init_displacement,init_velocity)

        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)

            self._solve(mp)
            current_analytical_displacement_y = A * cos(omega_D*time+theta) * exp(-delta*time)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0),current_analytical_displacement_y,delta=1e-3)

    def test_nodal_damping_explicit(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("sdof")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        mp.CreateNewProperties(1)

        self._add_variables(mp,explicit_dynamics=True)

        # Create node
        node = mp.CreateNewNode(1,0.0,0.0,0.0)
        _add_dofs(node)

        # Add bcs and initial values
        init_displacement = 0.1
        init_velocity = 0.0
        _set_dirichlet_bc(node,init_displacement,init_velocity)

        # Create element
        mass,stiffness,damping = _set_material_properties()
        _create_element(mp,mass,stiffness,damping)

        # Time integration parameters
        dt = 0.005
        time = 0.0
        end_time = 5.0
        step = 0

        self._set_and_fill_buffer(mp,2,dt)

        # Parameters for analytical solution
        omega_D,delta,theta,A = _return_parameters_analytical_solution(stiffness,mass,damping,init_displacement,init_velocity)

        strategy_expl = _create_dynamic_explicit_strategy(mp,'central_differences')
        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)

            strategy_expl.Solve()
            current_analytical_displacement_y = A * cos(omega_D*time+theta) * exp(-delta*time)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0),current_analytical_displacement_y,delta=1e-3)

def _create_dynamic_explicit_strategy(mp,scheme_name):
        if (scheme_name=='central_differences'):
            scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(0.00,0.00,0.00)
        elif scheme_name=='multi_stage':
            scheme = StructuralMechanicsApplication.ExplicitMultiStageKimScheme(0.33333333333333333)

        strategy = StructuralMechanicsApplication.MechanicalExplicitStrategy(mp,scheme,0,0,1)
        strategy.SetEchoLevel(0)
        return strategy

def _add_dofs(node):
    node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
    node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
    node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)

def _set_dirichlet_bc(node,init_displacement,init_velocity):
    node.Fix(KratosMultiphysics.DISPLACEMENT_X)
    node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
    node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,init_displacement)

def _create_element(mp,mass,stiffness,damping):
    prop = mp.GetProperties()[1]
    element = mp.CreateNewElement("NodalConcentratedDampedElement3D1N", 1, [1], prop)
    element.SetValue(KratosMultiphysics.NODAL_MASS,mass)
    element.SetValue(StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS,[0,stiffness,0])
    element.SetValue(StructuralMechanicsApplication.NODAL_DAMPING_RATIO,[0,damping,0])
    element.Initialize(mp.ProcessInfo)

def _return_parameters_analytical_solution(stiffness,mass,damping,init_displacement,init_velocity):
    omega = sqrt(stiffness/mass)
    D = damping / (2*mass*omega)
    omega_D = omega * sqrt(1-D*D)
    delta = damping / (2*mass)
    theta = atan(-(init_velocity+init_displacement*delta) / (omega_D*init_displacement))
    A = sqrt(init_displacement*init_displacement + ((init_velocity+init_displacement*delta) / omega_D)**2)
    return omega_D,delta,theta,A

def _set_material_properties():
    mass = 1.0
    stiffness = 10.0
    damping = 1.0
    return mass,stiffness,damping

if __name__ == '__main__':
    KratosUnittest.main()
