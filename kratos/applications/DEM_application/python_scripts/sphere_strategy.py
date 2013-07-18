import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

def Var_Translator(variable):

    if (variable == "OFF" or variable == "0"):
        variable = 0
    else:
        variable = 1

    return variable

def AddVariables(model_part, Param):

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
    model_part.AddNodalSolutionStepVariable(SQRT_OF_MASS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)
    model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS)
    model_part.AddNodalSolutionStepVariable(POISSON_RATIO)
    model_part.AddNodalSolutionStepVariable(LN_OF_RESTITUTION_COEFF)
    model_part.AddNodalSolutionStepVariable(PARTICLE_COHESION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_FRICTION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_TENSION)

    # ROTATION RELATED PROPERTIES

    if (Var_Translator(Param.RotationOption)):
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

    if (Var_Translator(Param.LimitSurfaceOption)):
        model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES)
        model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_ROTATE_SPRING_MOMENT)

    # FLAGS
    model_part.AddNodalSolutionStepVariable(GROUP_ID)            # Differencied groups for plotting, etc..
    model_part.AddNodalSolutionStepVariable(ERASE_FLAG)
    # ONLY VISUALITZATION

    if (Var_Translator(Param.PostGroupId)):
        model_part.AddNodalSolutionStepVariable(EXPORT_GROUP_ID)

    if (Var_Translator(Param.PostExportId)):
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

    def __init__(self, model_part, Param):

        # Initialization of member variables

        # SIMULATION FLAGS        
        self.virtual_mass_OPTION            = Var_Translator(Param.VirtualMassOption)
        self.critical_time_OPTION           = Var_Translator(Param.AutoReductionOfTimeStepOption)
        self.trihedron_OPTION               = Var_Translator(Param.TrihedronOption)
        self.rotation_OPTION                = Var_Translator(Param.RotationOption)
        self.bounding_box_OPTION            = Var_Translator(Param.BoundingBoxOption)
        self.fix_velocities                 = Var_Translator(Param.FixVelocitiesOption)
        self.limit_surface_OPTION           = Var_Translator(Param.LimitSurfaceOption)
        self.clean_init_indentation_OPTION  = Var_Translator(Param.CleanIndentationsOption)
        self.homogeneous_material_OPTION    = Var_Translator(Param.HomogeneousMaterialOption)
        self.global_variables_OPTION        = Var_Translator(Param.GlobalVariablesOption)
        self.Non_Linear_Option              = Var_Translator(Param.NonLinearNormalElasticOption)
        self.contact_mesh_OPTION            = Var_Translator(Param.ContactMeshOption)
        self.search_radius_extension        = Var_Translator(Param.SearchRadiusExtension)
        self.MoveMeshFlag                   = True
        self.deactivate_search              = 0
        self.case_OPTION                    = 3

        # MODEL
        self.model_part                      = model_part

        # BOUNDING_BOX
        self.enlargement_factor             = Param.BoundingBoxEnlargementFactor

        # BOUNDARY
        self.surface_normal_dir             = Vector(3)
        self.surface_normal_dir[0]          = Param.SurfaceNormalDirX
        self.surface_normal_dir[1]          = Param.SurfaceNormalDirY
        self.surface_normal_dir[2]          = Param.SurfaceNormalDirZ
        self.surface_point_coor             = Vector(3)
        self.surface_point_coor[0]          = Param.SurfacePointCoorX
        self.surface_point_coor[1]          = Param.SurfacePointCoorY
        self.surface_point_coor[2]          = Param.SurfacePointCoorZ
        self.surface_friction_angle         = Param.SurfaceFrictionAngle
        self.cylinder_radius                = Param.CylinderRadius

        # GLOBAL PHYSICAL ASPECTS
        self.gravity                        = Vector(3)
        self.gravity[0]                     = Param.GravityX
        self.gravity[1]                     = Param.GravityY
        self.gravity[2]                     = Param.GravityZ

        # GLOBAL MATERIAL PROPERTIES
        self.nodal_mass_coeff               = Param.VirtualMassCoefficient
        self.magic_factor                   = Param.MagicFactor

        if (self.global_variables_OPTION):
            self.global_kn                  = Param.GlobalKn
            self.global_kt                  = Param.GlobalKt

        if (Param.NormalForceCalculationType == "Linear"):
            self.force_calculation_type_id  = 0

        elif (Param.NormalForceCalculationType == "Hertz"):
            self.force_calculation_type_id  = 1

        else:

            raise 'Specified NormalForceCalculationType is not defined'

        if (self.Non_Linear_Option):
            self.C1                         = Param.C1
            self.C2                         = Param.C2
            self.N1                         = Param.N1
            self.N2                         = Param.N2

        if (Param.NormalDampingType == "ViscDamp"):

            if (Param.TangentialDampingType == "ViscDamp"):
                self.damp_id                = 11

            else:
                self.damp_id                = 10
        else:

            if (Param.TangentialDampingType == "ViscDamp"):
                self.damp_id                = 1

            else:
                self.damp_id                = 0

        if (Param.RotaDampingType == "LocalDamp"):
            self.rota_damp_id               = 1

        elif (Param.RotaDampingType == "RollingFric"):
            self.rota_damp_id               = 2

        else:
            self.rota_damp_id               = 0

        self.tau_zero                       = Param.TauZero
        self.sigma_max                      = Param.SigmaMax
        self.sigma_min                      = Param.SigmaMin

        # PRINTING VARIABLES
        self.print_export_id                = Var_Translator(Param.PostExportId)
        self.print_group_id                 = Var_Translator(Param.PostGroupId)
        self.print_radial_displacement      = Var_Translator(Param.PostRadialDisplacement)

        # TIME RELATED PARAMETERS
        self.delta_time                     = Param.MaxTimeStep
        self.max_delta_time                 = Param.MaxTimeStep
        self.final_time                     = Param.FinalTime

        # RESOLUTION METHODS AND PARAMETERS

        if (Param.TimeStepsPerSearchStep < 1):

            raise 'Variable TimeStepsPerSearchStep must be an integer, grater or equal to 1. The current input value is ', Param.TimeStepsPerSearchStep

        elif (not isinstance(Param.TimeStepsPerSearchStep, int)):

            print 'Variable TimeStepsPerSearchStep is not an integer. Its input value is ', Param.TimeStepsPerSearchStep, 'Rounding up to ', int(Param.TimeStepsPerSearchStep)

            self.n_step_search              = int(Param.TimeStepsPerSearchStep)

        else:
            self.n_step_search              = int(Param.TimeStepsPerSearchStep)

        if (self.deactivate_search):
            self.n_step_search              = sys.maxint

        self.safety_factor                  = Param.DeltaTimeSafetyFactor # For critical time step
        self.create_and_destroy             = ParticleCreatorDestructor()

        # STRATEGIES
        self.search_strategy                = OMP_DEMSearch()

        if (Param.IntegrationScheme == 'forward_euler'):
            self.time_scheme                = ForwardEulerScheme()

        elif (Param.IntegrationScheme == 'mid_point_rule'):
            self.time_scheme                = MidPointScheme()

        elif (Param.IntegrationScheme == 'const_average_acc'):
            self.time_scheme                = ConstAverageAccelerationScheme()

        else:

            print('Specified IntegrationScheme is not defined')

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
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_6, self.fix_velocities)
        self.model_part.ProcessInfo.SetValue(GLOBAL_VARIABLES_OPTION, self.global_variables_OPTION)
        self.model_part.ProcessInfo.SetValue(UNIFORM_MATERIAL_OPTION, self.homogeneous_material_OPTION)
        self.model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0);
        self.model_part.ProcessInfo.SetValue(TOTAL_CONTACTS, 0);
        self.model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_OPTION);

        # TOLERANCES
        self.model_part.ProcessInfo.SetValue(DISTANCE_TOLERANCE, 0);

        # BOUNDARY
        self.model_part.ProcessInfo.SetValue(LIMIT_SURFACE_OPTION, self.limit_surface_OPTION)
        self.model_part.ProcessInfo.SetValue(SURFACE_NORMAL_DIR, self.surface_normal_dir)
        self.model_part.ProcessInfo.SetValue(SURFACE_POINT_COOR, self.surface_point_coor)
        self.model_part.ProcessInfo.SetValue(SURFACE_FRICC, self.surface_friction_angle)
        self.model_part.ProcessInfo.SetValue(CYLINDER_RADIUS, self.cylinder_radius)

        # GLOBAL PHISICAL ASPECTS
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)
        self.model_part.ProcessInfo.SetValue(DEM_MAGIC_FACTOR, self.magic_factor)

        # GLOBAL MATERIAL PROPERTIES

        self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)
     
        if (self.global_variables_OPTION):
            self.model_part.ProcessInfo.SetValue(GLOBAL_KN, self.global_kn)
            self.model_part.ProcessInfo.SetValue(GLOBAL_KT, self.global_kt)

        # SEARCH-RELATED
        self.model_part.ProcessInfo.SetValue(SEARCH_RADIUS_EXTENSION, self.search_radius_extension)

        # PRINTING VARIABLES
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_8, self.print_group_id) # Reserved for: Export Print Group ID
        self.model_part.ProcessInfo.SetValue(FORCE_CALCULATION_TYPE, self.force_calculation_type_id)
        self.model_part.ProcessInfo.SetValue(DAMP_TYPE, self.damp_id)
        self.model_part.ProcessInfo.SetValue(ROTA_DAMP_TYPE, self.rota_damp_id)

        # TIME RELATED PARAMETERS
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)
        self.model_part.ProcessInfo.SetValue(FINAL_SIMULATION_TIME, self.final_time)
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_7, 0) # int(self.time_step_percentage_fix_velocities * (self.final_time / self.delta_time))) # Reserved for timestep fix_velocities

        # RESOLUTION METHODS AND PARAMETERS
        # Creating the solution strategy

        self.solver = ExplicitSolverStrategy(self.model_part, self.enlargement_factor, self.max_delta_time, self.n_step_search, self.safety_factor,
                                             self.MoveMeshFlag, self.time_scheme, self.search_strategy)

        self.solver.Initialize() # Calls the solver Initialized function (initializes all elements and performs other necessary tasks before iterating)

    #######################################################################

    def Solve(self):
        (self.solver).Solve()

    #######################################################################
