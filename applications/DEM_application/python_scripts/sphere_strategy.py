import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

def Var_Translator(variable):

    if (variable == "OFF" or variable == "0"):
        variable = 0
    else:
        variable = 1

    return variable

def AddVariables(ModelPart, Param):

    # KINEMATIC
    ModelPart.AddNodalSolutionStepVariable(DISPLACEMENT)
    ModelPart.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
    ModelPart.AddNodalSolutionStepVariable(VELOCITY)
    ModelPart.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
    ModelPart.AddNodalSolutionStepVariable(DELTA_ROTA_DISPLACEMENT)
    ModelPart.AddNodalSolutionStepVariable(ORIENTATION_REAL)
    ModelPart.AddNodalSolutionStepVariable(ORIENTATION_IMAG)
    ModelPart.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)

    # FORCES
    ModelPart.AddNodalSolutionStepVariable(RHS)
    ModelPart.AddNodalSolutionStepVariable(TOTAL_FORCES)
    ModelPart.AddNodalSolutionStepVariable(DAMP_FORCES)
    ModelPart.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
    ModelPart.AddNodalSolutionStepVariable(APPLIED_FORCE)

    # BASIC PARTICLE PROPERTIES
    ModelPart.AddNodalSolutionStepVariable(RADIUS)
    ModelPart.AddNodalSolutionStepVariable(NODAL_MASS)
    ModelPart.AddNodalSolutionStepVariable(SQRT_OF_MASS)
    ModelPart.AddNodalSolutionStepVariable(PARTICLE_DENSITY)
    ModelPart.AddNodalSolutionStepVariable(YOUNG_MODULUS)
    ModelPart.AddNodalSolutionStepVariable(POISSON_RATIO)
    ModelPart.AddNodalSolutionStepVariable(LN_OF_RESTITUTION_COEFF)
    ModelPart.AddNodalSolutionStepVariable(PARTICLE_COHESION)
    ModelPart.AddNodalSolutionStepVariable(PARTICLE_FRICTION)
    ModelPart.AddNodalSolutionStepVariable(PARTICLE_TENSION)

    # ROTATION RELATED PROPERTIES

    if (Var_Translator(Param.RotationOption)):
        ModelPart.AddNodalSolutionStepVariable(PARTICLE_INERTIA)
        ModelPart.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA)
        ModelPart.AddNodalSolutionStepVariable(PARTICLE_ROTATION_DAMP_RATIO)
        ModelPart.AddNodalSolutionStepVariable(ROLLING_FRICTION)

    # OTHER PROPERTIES
    ModelPart.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD
    ModelPart.AddNodalSolutionStepVariable(PARTICLE_CONTINUUM)  # Continuum group
    ModelPart.AddNodalSolutionStepVariable(REPRESENTATIVE_VOLUME)
    ModelPart.AddNodalSolutionStepVariable(MAX_INDENTATION)

    # LOCAL AXIS
    ModelPart.AddNodalSolutionStepVariable(EULER_ANGLES)

    # BOUNDARY SURFACE

    if (Var_Translator(Param.LimitSurfaceOption)):
        ModelPart.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES)
        ModelPart.AddNodalSolutionStepVariable(PARTICLE_SURFACE_ROTATE_SPRING_MOMENT)

    # FLAGS
    ModelPart.AddNodalSolutionStepVariable(GROUP_ID)            # Differencied groups for plotting, etc..

    # ONLY VISUALITZATION

    if (Var_Translator(Param.PostGroupId)):
        ModelPart.AddNodalSolutionStepVariable(EXPORT_GROUP_ID)

    if (Var_Translator(Param.PostExportId)):
        ModelPart.AddNodalSolutionStepVariable(EXPORT_ID)

    print "Variables for the explicit solver added correctly"

def AddDofs(ModelPart):

    for node in ModelPart.Nodes:
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

    def __init__(self, ModelPart, Param):

        # Initialization of member variables

        # SIMULATION FLAGS        
        self.virtual_mass_OPTION            = Var_Translator(Param.VirtualMassOption)
        self.critical_time_OPTION           = Var_Translator(Param.AutoReductionOfTimeStepOption)
        self.trihedron_OPTION               = Var_Translator(Param.TrihedronOption)
        self.rotation_OPTION                = Var_Translator(Param.RotationOption)
        self.rotation_spring_OPTION         = Var_Translator(Param.RotationalSpringOption)
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
        self.ModelPart                      = ModelPart

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
        self.create_and_destroy             = particle_destructor_and_constructor()

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
        self.ModelPart.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, self.virtual_mass_OPTION)
        self.ModelPart.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, self.critical_time_OPTION)
        self.ModelPart.ProcessInfo.SetValue(CASE_OPTION, self.case_OPTION)
        self.ModelPart.ProcessInfo.SetValue(TRIHEDRON_OPTION, self.trihedron_OPTION)
        self.ModelPart.ProcessInfo.SetValue(ROTATION_OPTION, self.rotation_OPTION)
        self.ModelPart.ProcessInfo.SetValue(BOUNDING_BOX_OPTION, self.bounding_box_OPTION)
        self.ModelPart.ProcessInfo.SetValue(INT_DUMMY_6, self.fix_velocities)
        self.ModelPart.ProcessInfo.SetValue(GLOBAL_VARIABLES_OPTION, self.global_variables_OPTION)
        self.ModelPart.ProcessInfo.SetValue(UNIFORM_MATERIAL_OPTION, self.homogeneous_material_OPTION)
        self.ModelPart.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0);
        self.ModelPart.ProcessInfo.SetValue(TOTAL_CONTACTS, 0);
        self.ModelPart.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_OPTION);
        self.ModelPart.ProcessInfo.SetValue(ROTATION_SPRING_OPTION, self.rotation_spring_OPTION);

        # TOLERANCES
        self.ModelPart.ProcessInfo.SetValue(DISTANCE_TOLERANCE, 0);

        # BOUNDARY
        self.ModelPart.ProcessInfo.SetValue(LIMIT_SURFACE_OPTION, self.limit_surface_OPTION)
        self.ModelPart.ProcessInfo.SetValue(SURFACE_NORMAL_DIR, self.surface_normal_dir)
        self.ModelPart.ProcessInfo.SetValue(SURFACE_POINT_COOR, self.surface_point_coor)
        self.ModelPart.ProcessInfo.SetValue(SURFACE_FRICC, self.surface_friction_angle)

        # GLOBAL PHISICAL ASPECTS
        self.ModelPart.ProcessInfo.SetValue(GRAVITY, self.gravity)
        self.ModelPart.ProcessInfo.SetValue(DEM_MAGIC_FACTOR, self.magic_factor)

        # GLOBAL MATERIAL PROPERTIES

        self.ModelPart.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)
     
        if (self.global_variables_OPTION):
            self.ModelPart.ProcessInfo.SetValue(GLOBAL_KN, self.global_kn)
            self.ModelPart.ProcessInfo.SetValue(GLOBAL_KT, self.global_kt)

        # SEARCH-RELATED
        self.ModelPart.ProcessInfo.SetValue(SEARCH_RADIUS_EXTENSION, self.search_radius_extension)

        # PRINTING VARIABLES
        self.ModelPart.ProcessInfo.SetValue(INT_DUMMY_8, self.print_group_id) # Reserved for: Export Print Group ID
        self.ModelPart.ProcessInfo.SetValue(FORCE_CALCULATION_TYPE, self.force_calculation_type_id)
        self.ModelPart.ProcessInfo.SetValue(DAMP_TYPE, self.damp_id)
        self.ModelPart.ProcessInfo.SetValue(ROTA_DAMP_TYPE, self.rota_damp_id)

        # TIME RELATED PARAMETERS
        self.ModelPart.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)
        self.ModelPart.ProcessInfo.SetValue(FINAL_SIMULATION_TIME, self.final_time)
        self.ModelPart.ProcessInfo.SetValue(INT_DUMMY_7, 0) # int(self.time_step_percentage_fix_velocities * (self.final_time / self.delta_time))) # Reserved for timestep fix_velocities

        # RESOLUTION METHODS AND PARAMETERS
        # Creating the solution strategy

        self.solver = ExplicitSolverStrategy(self.ModelPart, self.enlargement_factor, self.max_delta_time, self.n_step_search, self.safety_factor,
                                             self.MoveMeshFlag, self.time_scheme, self.search_strategy)

        self.solver.Initialize() # Calls the solver Initialized function (initializes all elements and performs other necessary tasks before iterating)

    #######################################################################

    def Solve(self):
        (self.solver).Solve()

    #######################################################################
