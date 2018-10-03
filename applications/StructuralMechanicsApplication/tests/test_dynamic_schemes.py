from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, cos, sin

class DynamicSchemesTests(KratosUnittest.TestCase):
    def setUp(self):
        pass
    

    def _base_spring_test_pseudo_static_scheme(self, scheme_name = "pseudo_static", buffer_size = 2, dt = 5.0e-3, beta = 0):
        mp = KratosMultiphysics.ModelPart("sdof")
        add_variables(mp, scheme_name)

        # Setting beta
        mp.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_BETA] = beta

        # Create node
        node = mp.CreateNewNode(1,0.0,0.0,0.0)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)

        #add bcs and initial values
        init_velocity = 0.1
        node.Fix(KratosMultiphysics.DISPLACEMENT_X)
        node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        node.Fix(KratosMultiphysics.VELOCITY_Y)
        node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0,init_velocity)

        #create element
        element = mp.CreateNewElement("NodalConcentratedElement3D1N", 1, [1], mp.GetProperties()[1])
        mass = 1.0
        stiffness = 0.0
        element.SetValue(KratosMultiphysics.NODAL_MASS, mass)
        element.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS, [0, stiffness,0])

        #time integration parameters
        time = 0.0
        end_time = 0.05
        step = 0

        set_and_fill_buffer(mp,buffer_size,dt)

        self.strategy = create_solver(mp, scheme_name)

        current_disp = 0.0
        current_vel = init_velocity
        # Solve the problem
        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)
            mp.ProcessInfo[KratosMultiphysics.STEP] = step

            self.strategy.Solve()

            current_analytical_displacement_y = current_disp + current_vel * dt - current_vel * beta * mass * dt/2.0
            current_disp = current_analytical_displacement_y
            current_analytical_velocity_y = current_vel - current_vel * beta * mass
            current_vel = current_analytical_velocity_y

            ## ASSERT
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0), current_analytical_displacement_y, delta=1e-3)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0), current_analytical_velocity_y, delta=1e-3)

    def _base_spring_test_dynamic_schemes(self, scheme_name = "bossak", buffer_size = 2, dt = 5.0e-3):
        mp = KratosMultiphysics.ModelPart("sdof")
        add_variables(mp, scheme_name)

        # Create node
        node = mp.CreateNewNode(1,0.0,0.0,0.0)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)

        #add bcs and initial values
        init_displacement = 0.1
        node.Fix(KratosMultiphysics.DISPLACEMENT_X)
        node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,init_displacement)

        #create element
        element = mp.CreateNewElement("NodalConcentratedElement3D1N", 1, [1], mp.GetProperties()[1])
        mass = 1.0
        stiffness = 10.0
        element.SetValue(KratosMultiphysics.NODAL_MASS, mass)
        element.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS, [0, stiffness,0])

        #time integration parameters
        time = 0.0
        end_time = 0.05
        step = 0

        set_and_fill_buffer(mp,buffer_size,dt)

        #parameters for analytical solution
        omega = sqrt(stiffness/mass)
        A = init_displacement

        self.strategy = create_solver(mp, scheme_name)

        # Fill buffer solution
        for i in range(buffer_size):
            time = time + dt
            step = i + 1
            mp.CloneTimeStep(time)
            mp.ProcessInfo[KratosMultiphysics.STEP] = step
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, A * cos(omega*time))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0, - A * omega * sin(omega*time))
            node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION_Y, 0, - A * omega**2 * cos(omega*time))

        # Solve the problem
        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)
            mp.ProcessInfo[KratosMultiphysics.STEP] = step

            self.strategy.Solve()

            current_analytical_displacement_y = A * cos(omega*time)
            current_analytical_velocity_y = - A *omega * sin(omega*time)
            current_analytical_acceleration_y = - A * omega * omega * cos(omega*time)

            # ASSERT
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0), current_analytical_displacement_y, delta=1e-3)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0), current_analytical_velocity_y, delta=1e-3)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION_Y,0), current_analytical_acceleration_y, delta=1e-3)

    def _base_fall_test_dynamic_schemes(self, scheme_name = "bossak", buffer_size = 2, dt = 1.0e-2):
        mp = KratosMultiphysics.ModelPart("sdof")
        add_variables(mp, scheme_name)

        if (scheme_name == "explicit"):
            # Create node
            node = mp.CreateNewNode(1,0.0,0.0,0.0)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)
            second_node = mp.CreateNewNode(2,1.0,0.0,0.0)
            second_node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
            second_node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
            second_node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)

            #add bcs and initial values
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,0.0)
            second_node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            second_node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
            second_node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,0.0)

            #create element
            prop = mp.GetProperties()[1]
            prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, StructuralMechanicsApplication.TrussConstitutiveLaw())
            prop.SetValue(KratosMultiphysics.DENSITY, 1.0)
            prop.SetValue(StructuralMechanicsApplication.CROSS_AREA, 1.0)
            prop.SetValue(KratosMultiphysics.YOUNG_MODULUS, 1.0)
            prop.SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2, 0.0)
            element = mp.CreateNewElement("TrussElement3D2N", 1, [1, 2], prop)
            gravity = -9.81
            node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_Y,0,gravity)
            second_node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_Y,0,gravity)
        else:
            # Create node
            node = mp.CreateNewNode(1,0.0,0.0,0.0)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)

            #add bcs and initial values
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,0.0)

            #create element
            element = mp.CreateNewElement("NodalConcentratedElement3D1N", 1, [1], mp.GetProperties()[1])
            mass = 1.0
            stiffness = 0.0
            element.SetValue(KratosMultiphysics.NODAL_MASS, mass)
            element.SetValue(StructuralMechanicsApplication.NODAL_STIFFNESS, [0, stiffness,0])
            gravity = -9.81
            node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_Y,0,gravity)

        #time integration parameters
        time = 0.0
        end_time = 1.0e-1
        step = 0

        set_and_fill_buffer(mp,buffer_size,dt)

        #parameters for analytical solution
        self.strategy = create_solver(mp, scheme_name)

        current_analytical_displacement_y = 0.0
        current_analytical_velocity_y = 0.0
        current_analytical_acceleration_y = 0.0

        # Fill buffer solution
        for i in range(buffer_size):
            time = time + dt
            step = i + 1
            mp.CloneTimeStep(time)
            mp.ProcessInfo[KratosMultiphysics.STEP] = step

            current_analytical_displacement_y += current_analytical_velocity_y * dt + 0.25 * (current_analytical_acceleration_y + gravity) * dt**2
            current_analytical_velocity_y += dt * 0.5 * (current_analytical_acceleration_y + gravity)
            current_analytical_acceleration_y = gravity

            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, current_analytical_displacement_y)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0, current_analytical_velocity_y)
            node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION_Y, 0, current_analytical_acceleration_y)

        # Solve the problem
        while(time <= end_time):
            time = time + dt
            step = step + 1
            mp.CloneTimeStep(time)
            mp.ProcessInfo[KratosMultiphysics.STEP] = step

            self.strategy.Solve()
            current_analytical_displacement_y += current_analytical_velocity_y * dt + 0.25 * (current_analytical_acceleration_y + gravity) * dt**2
            current_analytical_velocity_y += dt * 0.5 * (current_analytical_acceleration_y + gravity)
            current_analytical_acceleration_y = gravity

            # ASSERT
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION_Y,0), current_analytical_acceleration_y, delta=1e-3)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0), current_analytical_velocity_y, delta=1e-3)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0), current_analytical_displacement_y, delta=1e-3)

    def test_spring_bossak_scheme(self):
        current_model = KratosMultiphysics.Model()
        self._base_spring_test_dynamic_schemes(current_model,"bossak", 2)

    def test_spring_newmark_scheme(self):
        current_model = KratosMultiphysics.Model()
        self._base_spring_test_dynamic_schemes(current_model,"newmark", 2)

    def test_spring_backward_euler_scheme(self):
        current_model = KratosMultiphysics.Model()
        self._base_spring_test_dynamic_schemes(current_model,"backward_euler", 2, 4.0e-3)

    def test_spring_bdf2_scheme(self):
        current_model = KratosMultiphysics.Model()
        self._base_spring_test_dynamic_schemes(current_model,"bdf2", 3)

    def test_fall_bossak_scheme(self):
        current_model = KratosMultiphysics.Model()
        self._base_fall_test_dynamic_schemes(current_model,"bossak", 2)

    def test_fall_newmark_scheme(self):
        current_model = KratosMultiphysics.Model()
        self._base_fall_test_dynamic_schemes(current_model,"newmark", 2)

    def test_fall_backward_euler_scheme(self):
        current_model = KratosMultiphysics.Model()
        self._base_fall_test_dynamic_schemes(current_model,"backward_euler", 2, 2.0e-3)

    def test_fall_bdf2_scheme(self):
        current_model = KratosMultiphysics.Model()
        self._base_fall_test_dynamic_schemes(current_model,"bdf2", 3)

    def test_fall_explicit_scheme(self):
        self._base_fall_test_dynamic_schemes("explicit", 2, 1.0e-5)

    def test_spring_test_pseudo_static_scheme(self):
        self._base_spring_test_pseudo_static_scheme("pseudo_static", 2, 1.0e-2)

    def test_spring_test_pseudo_static_with_damping_scheme(self):
        self._base_spring_test_pseudo_static_scheme("pseudo_static", 2, 1.0e-2, 1.0)

def set_and_fill_buffer(mp,buffer_size,delta_time):
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

def add_variables(mp, scheme_name = "bossak"):
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
    if (scheme_name == "explicit"):
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)


def create_solver(mp, scheme_name):
    # Define a minimal newton raphson dynamic solver
    linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
    builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
    if (scheme_name == "newmark"):
        scheme = KratosMultiphysics.ResidualBasedNewmarkDisplacementScheme()
    elif (scheme_name == "pseudo_static"):
        scheme = KratosMultiphysics.ResidualBasedPseudoStaticDisplacementScheme(StructuralMechanicsApplication.RAYLEIGH_BETA)
    elif (scheme_name == "backward_euler"):
        scheme = KratosMultiphysics.ResidualBasedBDFDisplacementScheme(1)
    elif (scheme_name == "bdf2"):
        scheme = KratosMultiphysics.ResidualBasedBDFDisplacementScheme(2)
    elif (scheme_name == "explicit"):
        dynamic_settings = KratosMultiphysics.Parameters("""
        {
            "time_step_prediction_level": 0,
            "max_delta_time": 1.0e-5,
            "fraction_delta_time": 0.9
        }
        """)
        scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(dynamic_settings)
    else:
        damp_factor_m = 0.0
        scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m)


    # Creating solver
    compute_reactions = False
    reform_step_dofs = True
    move_mesh_flag = True

    # Explicit solver
    if (scheme_name == "explicit"):
        strategy = StructuralMechanicsApplication.MechanicalExplicitStrategy(mp,
                                                                        scheme,
                                                                        compute_reactions,
                                                                        reform_step_dofs,
                                                                        move_mesh_flag)
    else: # Implicit solver otherwise
        max_iters = 10
        # Convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-14,1e-20)
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-4,1e-9)
        convergence_criterion.SetEchoLevel(0)
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

if __name__ == '__main__':
    KratosUnittest.main()
