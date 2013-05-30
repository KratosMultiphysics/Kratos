
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
    #model_part.AddNodalSolutionStepVariable(APPLIED_FORCE)

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

    # LOCAL AXIS
    model_part.AddNodalSolutionStepVariable(EULER_ANGLES)

    #SURFACE
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

class ExplicitStrategy:

    def __init__(self, model_part, domain_size):

        # Initialization of member variables

        # MODEL
        self.model_part                     = model_part
        self.domain_size                    = domain_size

        # BOUNDARY
        self.limit_surface_OPTION           = 0 #its 1/0 xapuza
        self.surface_normal_dir             = Vector(3)
        self.surface_normal_dir[0]          = 0.0
        self.surface_normal_dir[1]          = 1.0
        self.surface_normal_dir[2]          = 0.0
        self.surface_point_coor             = Vector(3)
        self.surface_point_coor[0]          = 0.0
        self.surface_point_coor[1]          = 0.0
        self.surface_point_coor[2]          = 0.0
        self.surface_friction_angle         = 45

        # GLOBAL PHISICAL ASPECTS
        self.damping_ratio                  = 0.00
        self.penalty_factor                 = 10.00
        self.gravity                        = Vector(3)
        self.gravity[0]                     = 0.0
        self.gravity[1]                     = - 9.81
        self.gravity[2]                     = 0.0
        self.nodal_mass_coeff               = 0.0
        self.magic_factor                   = 1.0
        self.global_kn                      = 1000.0
        self.global_kt                      = 1000.0

        # SIMULATION FLAGS
        self.MoveMeshFlag                   = True
        self.virtual_mass_OPTION            = 0  #its 1/0 xapuza
        self.critical_time_OPTION           = 0  #its 1/0 xapuza
        self.case_OPTION                    = 3  #aixo es una xapuza fins que pooyan permeti bools a pyton o tinguis flags.
        self.trihedron_OPTION               = 0
        self.rotation_OPTION                = 0  #its 1/0 xapuza
        self.rotation_spring_OPTION         = 0  #its 1/0 xapuza
        self.bounding_box_OPTION            = 0  #its 1/0 xapuza
        self.activate_search                = 1  #its 1/0 xapuza
        self.fix_velocities                 = 0
        self.global_variables_OPTION        = 0  #its 1/0 xapuza

        # PRINTING VARIABLES
        self.print_export_id                = 0
        self.print_export_skin_sphere       = 0
        self.force_calculation_type_id      = 1
        self.damp_id                        = 1
        self.rota_damp_id                   = 0
        self.dummy_switch                   = 0

        # TIME RELATED PARAMETERS
        self.delta_time                     = 0.00001
        self.max_delta_time                 = 0.05
        self.fraction_delta_time            = 0.90
        self.final_time                     = 3.0
        self.time_increasing_ratio          = 15 # Percentage%

        # RESOLUTION METHODS AND PARAMETERS
        self.enlargement_factor             = 1
        self.n_step_search                  = 1
        self.safety_factor                  = 1.0 # For critical time step
        self.create_and_destroy             = particle_destructor_and_constructor();

        # STRATEGIES
    self.time_scheme                    = ForwardEulerScheme();
    self.search_strategy                = OMP_DEMSearch();

    ######################################################################

    def Initialize(self):

        # Setting ProcessInfo variables

        # BOUNDARY
        self.model_part.ProcessInfo.SetValue(LIMIT_SURFACE_OPTION, self.limit_surface_OPTION)
        self.model_part.ProcessInfo.SetValue(SURFACE_NORMAL_DIR, self.surface_normal_dir)
        self.model_part.ProcessInfo.SetValue(SURFACE_POINT_COOR, self.surface_point_coor)
        self.model_part.ProcessInfo.SetValue(SURFACE_FRICC, self.surface_friction_angle)

        # GLOBAL PHISICAL ASPECTS
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)
        self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)
        self.model_part.ProcessInfo.SetValue(DEM_MAGIC_FACTOR, self.magic_factor)
        self.model_part.ProcessInfo.SetValue(GLOBAL_KN, self.global_kn)
        self.model_part.ProcessInfo.SetValue(GLOBAL_KT, self.global_kt)

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
        self.model_part.ProcessInfo.SetValue(UNIFORM_MATERIAL_OPTION, 1)
        self.model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0);
        self.model_part.ProcessInfo.SetValue(TOTAL_CONTACTS, 0);

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

        self.solver = ExplicitSolverStrategy(self.model_part, self.enlargement_factor, self.fraction_delta_time, self.n_step_search, self.safety_factor,
                                             self.MoveMeshFlag, self.time_scheme, self.search_strategy)

        self.solver.Initialize() # Calls the solver Initialized function (initializes all elements and performs other necessary tasks before iterating)

    #######################################################################

    def Initial_Critical_Time(self):
        (self.solver).InitialCriticalTime()

    #######################################################################

    def Solve(self):
        (self.solver).Solve()

    #######################################################################
