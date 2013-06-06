
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from DEM_explicit_solver_var import *

def AddVariables(model_part):

    # KINEMATIC
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
    model_part.AddNodalSolutionStepVariable(DELTA_ROTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(ORIENTATION_REAL)
    model_part.AddNodalSolutionStepVariable(ORIENTATION_IMAG)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)

    # FORCES
    model_part.AddNodalSolutionStepVariable(RHS)
    model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
    model_part.AddNodalSolutionStepVariable(DAMP_FORCES)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
    model_part.AddNodalSolutionStepVariable(APPLIED_FORCE)

    # BASIC PARTICLE PROPERTIES
    model_part.AddNodalSolutionStepVariable(RADIUS)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)
    model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS)
    model_part.AddNodalSolutionStepVariable(POISSON_RATIO)
    model_part.AddNodalSolutionStepVariable(RESTITUTION_COEFF)
    model_part.AddNodalSolutionStepVariable(PARTICLE_COHESION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_FRICTION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_TENSION)

    # ROTATION RELATED PROPERTIES
    model_part.AddNodalSolutionStepVariable(PARTICLE_INERTIA)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA)
    model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_DAMP_RATIO)
    model_part.AddNodalSolutionStepVariable(ROLLING_FRICTION)

    # OTHER PROPERTIES
    model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD
    model_part.AddNodalSolutionStepVariable(PARTICLE_CONTINUUM)  # Continuum group
    model_part.AddNodalSolutionStepVariable(REPRESENTATIVE_VOLUME)
    model_part.AddNodalSolutionStepVariable(MAX_INDENTATION)

    # LOCAL AXIS
    model_part.AddNodalSolutionStepVariable(EULER_ANGLES)

    # BOUNDARY SURFACE
    model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES)
    model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_ROTATE_SPRING_MOMENT)

    # FLAGS
    model_part.AddNodalSolutionStepVariable(GROUP_ID)            # Differencied groups for plotting, etc..

    # ONLY VISUALITZATION

    if (print_export_id == "1"):
      model_part.AddNodalSolutionStepVariable(EXPORT_ID)

    print "Variables for the explicit solver added correctly"

def AddDofs(model_part):

    for node in model_part.Nodes:
        node.AddDof(DISPLACEMENT_X, REACTION_X);
        node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        node.AddDof(VELOCITY_X, REACTION_X);
        node.AddDof(VELOCITY_Y, REACTION_Y);
        node.AddDof(VELOCITY_Z, REACTION_Z);
        node.AddDof(ANGULAR_VELOCITY_X, REACTION_X);
        node.AddDof(ANGULAR_VELOCITY_Y, REACTION_Y);
        node.AddDof(ANGULAR_VELOCITY_Z, REACTION_Z);

    print "DOFs for the DEM solution added correctly"

def Var_Translator(variable):

    if (variable == "OFF" or variable == "0"):
        variable = 0
    else:
        variable = 1

    return variable

class ExplicitStrategy:

    def __init__(self, model_part, domain_size):

        # Initialization of member variables

        # SIMULATION FLAGS        
        self.virtual_mass_OPTION            = Var_Translator(VirtualMassOption)  #its 1/0 xapuza
        self.critical_time_OPTION           = Var_Translator(AutoReductionOfTimeStepOption)  #its 1/0 xapuza
        self.case_OPTION                    = 3  #aixo es una xapuza fins que pooyan permeti bools a pyton o tinguis flags.
        self.trihedron_OPTION               = Var_Translator(TrihedronOption)
        self.rotation_OPTION                = Var_Translator(RotationOption)
        self.rotation_spring_OPTION         = Var_Translator(RotationalSpringOption)  #its 1/0 xapuza
        self.bounding_box_OPTION            = Var_Translator(BoundingBoxOption)  #its 1/0 xapuza
        self.activate_search                = 1  #its 1/0 xapuza
        self.fix_velocities                 = Var_Translator(FixVelocities)
        self.limit_surface_OPTION           = Var_Translator(LimitSurfaceOption)  #its 1/0 xapuza
        self.clean_init_indentation_OPTION  = Var_Translator(CleanIndentationsOption)
        self.homogeneous_material_OPTION    = Var_Translator(HomogeneousMaterialOption)
        self.global_variables_OPTION        = Var_Translator(GlobalVariablesOption)
        self.Non_Linear_Option              = Var_Translator(NonLinearOption)
        self.stress_strain_operations       = Var_Translator(StressStrainOperations)
        self.contact_mesh_OPTION            = Var_Translator(ContactMeshOption)
        self.concrete_test_OPTION           = Var_Translator(ConcreteTestOption)
        self.concrete_test_OPTION           = Var_Translator(ConcreteTestOption)
        self.search_radius_extension        = Var_Translator(DeltaOption)
        self.MoveMeshFlag                   = True


        if (self.search_radius_extension and Var_Translator(ContinuumOption)):
            self.amplified_continuum_search_radius_extension = Var_Translator(ContinuumOption)

        if (Var_Translator(DeltaOption)):
            self.delta_OPTION               = True

        if (Var_Translator(ContinuumOption)):
            self.continuum_simulating_OPTION = True

        # MODEL
        self.model_part                     = model_part
        self.domain_size                    = domain_size

        # BOUNDARY
        self.surface_normal_dir             = Vector(3)
        self.surface_normal_dir[0]          = surface_normal_dir_x
        self.surface_normal_dir[1]          = surface_normal_dir_y
        self.surface_normal_dir[2]          = surface_normal_dir_z
        self.surface_point_coor             = Vector(3)
        self.surface_point_coor[0]          = surface_point_coor_x
        self.surface_point_coor[1]          = surface_point_coor_y
        self.surface_point_coor[2]          = surface_point_coor_z
        self.surface_friction_angle         = surface_friction_angle


        # GLOBAL PHISICAL ASPECTS
        self.gravity                        = Vector(3)
        self.gravity[0]                     = gravity_x
        self.gravity[1]                     = gravity_y
        self.gravity[2]                     = gravity_z

        # GLOBAL MATERIAL PROPERTIES
        self.nodal_mass_coeff               = VirtualMassCoefficient
        self.magic_factor                   = Var_Translator(MagicFactor)

        if (self.global_variables_OPTION):
            self.global_kn                  = global_kn
            self.global_kt                  = global_kt
            solver.global_kr                = global_kr
            solver.global_rn                = global_rn
            solver.global_rt                = global_rt
            solver.global_rr                = global_rr
            solver.global_fri_ang           = global_fri_ang

        if (NormalForceCalculation == "Linear"):
            self.force_calculation_type_id  = 0
        elif (NormalForceCalculation == "Hertz"):
            self.force_calculation_type_id  = 1

        if (self.Non_Linear_Option):
            solver.C1                       = C1
            solver.C2                       = C2
            solver.N1                       = N1
            solver.N2                       = N2

        if(NormalDampId == "ViscDamp"):
            if (TangentialDampId == "ViscDamp"):
                self.damp_id=damp_id        = 11
            else:
                self.damp_id=damp_id        = 10
        else:
            if (TangentialDampId == "ViscDamp"):
                self.damp_id                = 1
            else:
                self.damp_id                = 0

        if (RotaDampId == "LocalDamp"):
            self.rota_damp_id               = 1
        elif (RotaDampId == "RollingFric"):
            self.rota_damp_id               = 2
        else:
            self.rota_damp_id               = 0

        if (FailureCriterionOption == "Mohr-Coulomb"):
            self.failure_criterion_OPTION   = 1
        elif (FailureCriterionOption == "Uncoupled"):
            self.failure_criterion_OPTION   = 2

        self.tau_zero                       = TauZero
        self.sigma_max                      = SigmaMax
        self.sigma_min                      = SigmaMin
        self.internal_fricc                 = InternalFricc

        # CONCRETE TEST

        if (Var_Translator(TriaxialOption)):
          self.initial_pressure_time        = InitialTime
          self.time_increasing_ratio        = IncreasingTemporaily

        # PRINTING VARIABLES
        self.print_export_id                = Var_Translator(print_export_id)
        self.print_export_skin_sphere       = Var_Translator(print_export_skin_sphere)
        self.damp_id                        = 1
        self.print_radial_displacement      = Var_Translator(print_radial_displacement)

        # TIME RELATED PARAMETERS
        self.delta_time                     = max_time_step
        self.max_delta_time                 = max_time_step
        self.final_time                     = final_time
        self.time_increasing_ratio          = int(IncreasingTemporaily) # Percentage (%)

        # RESOLUTION METHODS AND PARAMETERS
        self.n_step_search                  = int(search_step)
        self.safety_factor                  = dt_safety_factor # For critical time step
        self.create_and_destroy             = particle_destructor_and_constructor()

        # STRATEGIES
        self.search_strategy                = OMP_DEMSearch()

        if (Integration_Scheme == 'forward_euler'):
            time_scheme = ForwardEulerScheme()
        elif (Integration_Scheme == 'mid_point_rule'):
            time_scheme = MidPointScheme()
        elif (Integration_Scheme == 'const_average_acc'):
            time_scheme = ConstAverageAccelerationScheme()
        else:
            print('scheme not defined')

    ######################################################################

    def Initialize(self):

        # Setting ProcessInfo variables

        # SIMULATION FLAGS
        self.model_part.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, self.virtual_mass_OPTION)
        self.model_part.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, self.critical_time_OPTION)
        self.model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_OPTION)
        self.model_part.ProcessInfo.SetValue(TRIHEDRON_OPTION, self.trihedron_OPTION)
        self.model_part.ProcessInfo.SetValue(ROTATION_OPTION, self.rotation_OPTION)
        self.model_part.ProcessInfo.SetValue(BOUNDING_BOX_OPTION, self.bounding_box_OPTION)
        self.model_part.ProcessInfo.SetValue(ACTIVATE_SEARCH, self.activate_search)
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_6, self.fix_velocities) #reserved for fix_velocities
        self.model_part.ProcessInfo.SetValue(GLOBAL_VARIABLES_OPTION, self.global_variables_OPTION)
        self.model_part.ProcessInfo.SetValue(UNIFORM_MATERIAL_OPTION, self.homogeneous_material_OPTION)
        self.model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0);
        self.model_part.ProcessInfo.SetValue(TOTAL_CONTACTS, 0);
        self.model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_OPTION);
        self.model_part.ProcessInfo.SetValue(ROTATION_SPRING_OPTION, self.rotation_spring_OPTION);


        # TOLERANCES
        self.model_part.ProcessInfo.SetValue(DISTANCE_TOLERANCE, 0);

        # BOUNDARY
        self.model_part.ProcessInfo.SetValue(LIMIT_SURFACE_OPTION, self.limit_surface_OPTION)
        self.model_part.ProcessInfo.SetValue(SURFACE_NORMAL_DIR, self.surface_normal_dir)
        self.model_part.ProcessInfo.SetValue(SURFACE_POINT_COOR, self.surface_point_coor)
        self.model_part.ProcessInfo.SetValue(SURFACE_FRICC, self.surface_friction_angle)

        # GLOBAL PHISICAL ASPECTS
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)
        self.model_part.ProcessInfo.SetValue(DEM_MAGIC_FACTOR, self.magic_factor)

        # GLOBAL MATERIAL PROPERTIES

        if(self.homogeneous_material_OPTION == "ON"):
            self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)
            self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.magic_factor)
            self.model_part.ProcessInfo.SetValue(HISTORICAL_MIN_K, self.magic_factor)

        if (self.global_variables_OPTION == "ON"):
            self.model_part.ProcessInfo.SetValue(GLOBAL_KN, self.global_kn)
            self.model_part.ProcessInfo.SetValue(GLOBAL_KT, self.global_kt)

        # PRINTING VARIABLES
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_3, self.print_export_id) # Reserved for: Export Print Skin sphere
        self.model_part.ProcessInfo.SetValue(FORCE_CALCULATION_TYPE, self.force_calculation_type_id)
        self.model_part.ProcessInfo.SetValue(DAMP_TYPE, self.damp_id)
        self.model_part.ProcessInfo.SetValue(ROTA_DAMP_TYPE, self.rota_damp_id)
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_1, 0) # Reserved for: message when confinement ends.

        # TIME RELATED PARAMETERS
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)
        self.model_part.ProcessInfo.SetValue(FINAL_SIMULATION_TIME, self.final_time)
        self.model_part.ProcessInfo.SetValue(TIME_INCREASING_RATIO, self.time_increasing_ratio)
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_7, 0) # int(self.time_step_percentage_fix_velocities * (self.final_time / self.delta_time))) # Reserved for timestep fix_velocities

        # RESOLUTION METHODS AND PARAMETERS
        # Creating the solution strategy

        self.solver = ExplicitSolverStrategy(self.model_part, self.enlargement_factor, self.max_delta_time, self.n_step_search, self.safety_factor,
                                             self.MoveMeshFlag, self.time_scheme, self.search_strategy)

        self.solver.Initialize() # Calls the solver Initialized function (initializes all elements and performs other necessary tasks before iterating)

    #######################################################################

    def Initial_Critical_Time(self):
        (self.solver).InitialTimeStepCalculation()

    #######################################################################

    def Solve(self):
        (self.solver).Solve()

    #######################################################################
