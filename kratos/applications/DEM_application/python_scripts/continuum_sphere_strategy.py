from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *


def AddVariables(model_part, Param):

    # KINEMATIC
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
    model_part.AddNodalSolutionStepVariable(DELTA_ROTA_DISPLACEMENT)
#    model_part.AddNodalSolutionStepVariable(ORIENTATION_REAL)
#    model_part.AddNodalSolutionStepVariable(ORIENTATION_IMAG)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)

    # FORCES
    model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
    model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
    model_part.AddNodalSolutionStepVariable(DAMP_FORCES)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_FORCE)

    # BASIC PARTICLE PROPERTIES
    model_part.AddNodalSolutionStepVariable(RADIUS)
#    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(SQRT_OF_MASS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)
    model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS)
    model_part.AddNodalSolutionStepVariable(POISSON_RATIO)
    model_part.AddNodalSolutionStepVariable(LN_OF_RESTITUTION_COEFF)
    model_part.AddNodalSolutionStepVariable(PARTICLE_FRICTION)

    # ROTATION RELATED PROPERTIES
    if (Var_Translator(Param.RotationOption)):
        #model_part.AddNodalSolutionStepVariable(PARTICLE_INERTIA)
        model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA)
        model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_DAMP_RATIO)
        if( Var_Translator(Param.RollingFrictionOption)):
          model_part.AddNodalSolutionStepVariable(ROLLING_FRICTION)

    # OTHER PROPERTIES
    model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD
    model_part.AddNodalSolutionStepVariable(COHESIVE_GROUP)  # Continuum group
#    model_part.AddNodalSolutionStepVariable(REPRESENTATIVE_VOLUME)
    model_part.AddNodalSolutionStepVariable(MAX_INDENTATION)

    # LOCAL AXIS
    if (Param.PostEulerAngles == "1" or Param.PostEulerAngles == 1):
        model_part.AddNodalSolutionStepVariable(EULER_ANGLES)

# BOUNDARY SURFACE
#
#    if (Param.LimitSurfaceOption > 0):
#        model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES_1)
#    if (Param.LimitSurfaceOption > 1):
#        model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES_2)
#    if (Param.LimitSurfaceOption > 2):
#        model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES_3)
#    if (Param.LimitSurfaceOption > 3):
#        model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES_4)
#    if (Param.LimitSurfaceOption > 4):
#        model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES_5)
#
#    if (Param.LimitCylinderOption > 0):
#        model_part.AddNodalSolutionStepVariable(PARTICLE_CYLINDER_CONTACT_FORCES_1)
#    if (Param.LimitCylinderOption > 1):
#        model_part.AddNodalSolutionStepVariable(PARTICLE_CYLINDER_CONTACT_FORCES_2)
#    if (Param.LimitCylinderOption > 2):
#        model_part.AddNodalSolutionStepVariable(PARTICLE_CYLINDER_CONTACT_FORCES_3)
#    if (Param.LimitCylinderOption > 3):
#        model_part.AddNodalSolutionStepVariable(PARTICLE_CYLINDER_CONTACT_FORCES_4)
#    if (Param.LimitCylinderOption > 4):
#        model_part.AddNodalSolutionStepVariable(PARTICLE_CYLINDER_CONTACT_FORCES_5)
#
    # FLAGS
    model_part.AddNodalSolutionStepVariable(GROUP_ID)            # Differencied groups for plotting, etc..
#    model_part.AddNodalSolutionStepVariable(ERASE_FLAG)

    # OPTIMIZATION
    model_part.AddNodalSolutionStepVariable(VELOCITY_X_DOF_POS)
    model_part.AddNodalSolutionStepVariable(VELOCITY_Y_DOF_POS)
    model_part.AddNodalSolutionStepVariable(VELOCITY_Z_DOF_POS)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY_X_DOF_POS)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY_Y_DOF_POS)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY_Z_DOF_POS)
    model_part.AddNodalSolutionStepVariable(OLD_COORDINATES)

    # ONLY VISUALIZATION

    if (Var_Translator(Param.PostExportId)):
        model_part.AddNodalSolutionStepVariable(EXPORT_ID)
    if (Var_Translator(Param.PredefinedSkinOption) | Var_Translator(Param.ConcreteTestOption)):
        model_part.AddNodalSolutionStepVariable(EXPORT_SKIN_SPHERE)
        model_part.AddNodalSolutionStepVariable(PREDEFINED_SKIN)
    if (Var_Translator(Param.PostGroupId)):
        model_part.AddNodalSolutionStepVariable(EXPORT_GROUP_ID)

    print("Variables for the explicit solver added correctly")


def AddDofs(model_part):

    for node in model_part.Nodes:
        node.AddDof(DISPLACEMENT_X, REACTION_X)
        node.AddDof(DISPLACEMENT_Y, REACTION_Y)
        node.AddDof(DISPLACEMENT_Z, REACTION_Z)
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(ANGULAR_VELOCITY_X, REACTION_X);
        node.AddDof(ANGULAR_VELOCITY_Y, REACTION_Y);
        node.AddDof(ANGULAR_VELOCITY_Z, REACTION_Z);

    print("DOFs for the DEM solution added correctly")


def Var_Translator(variable):

    if (variable == "OFF" or variable == "0" or variable == 0):
        variable = 0
    else:
        variable = 1

    return variable


class ExplicitStrategy:

    def __init__(self, model_part, fem_model_part, creator_destructor, Param):

        # Initialization of member variables


        # SIMULATION FLAGS        
  
        self.critical_time_option           = Var_Translator(Param.AutoReductionOfTimeStepOption)   
        self.case_option                    = 3  
        self.trihedron_option               = Var_Translator(Param.PostEulerAngles)
        self.rotation_option                = Var_Translator(Param.RotationOption)
        self.bounding_box_option            = Var_Translator(Param.BoundingBoxOption)
        self.activate_search                = 1
        if(Var_Translator(Param.DontSearchUntilFailure)):
          print ("Search is not active until a bond is broken.")
          self.activate_search                = 0
          if (len(fem_model_part.Nodes)>0 or Param.ConcreteTestOption== "BTS" ):   #MSI. This activates the search since there are fem contact elements. however only the particle - fem search should be active.
            print ("WARNING: Search should be activated since there might contact with FEM.")

        self.fix_velocities_flag                 = 0       

        self.clean_init_indentation_option = Var_Translator(Param.CleanIndentationsOption)
        self.stress_strain_operations = Var_Translator(Param.StressStrainOperationsOption)
        self.MoveMeshFlag = True
        
        self.virtual_mass_option            = 0
        self.nodal_mass_coeff = Param.VirtualMassCoefficient
        if(self.nodal_mass_coeff != 1.00):
           self.virtual_mass_option            = 1
           
        
        self.delta_option = Var_Translator(Param.DeltaOption)

        self.continuum_simulating_option = Var_Translator(Param.ContinuumOption)
        self.dempack_option = Var_Translator(Param.Dempack)
        if(self.dempack_option):
            self.dempack_damping = Param.DempackDamping
            self.dempack_global_damping = Param.DempackGlobalDamping
            self.area_factor = Param.AreaFactor
        self.contact_mesh_option = Var_Translator(Var_Translator(Param.ContactMeshOption) & Var_Translator(Param.ContinuumOption))
        self.concrete_test_option = Param.ConcreteTestOption
        self.automatic_bounding_box_option = Var_Translator(Param.AutomaticBoundingBoxOption)

        self.search_tolerance = 0.0
        self.coordination_number = 10.0
        self.amplified_continuum_search_radius_extension = 1.0;

        if (Param.DeltaOption == "OFF"):
            self.delta_option = 0

        elif (Param.DeltaOption == "ABSOLUTE"):
            self.delta_option = 1
            self.search_tolerance = Param.SearchTolerance

        elif (Param.DeltaOption == "COORDINATION_NUMBER"):
            self.delta_option = 2
            self.coordination_number = Param.CoordinationNumber
            self.search_tolerance = 0.01 * Param.MeanRadius

        if (self.continuum_simulating_option):
            self.amplified_continuum_search_radius_extension = Param.AmplifiedSearchRadiusExtension;

        if(self.delta_option > 0):
            if(self.continuum_simulating_option):
                self.case_option = 2
            else:
                self.case_option = 1
        elif(not self.delta_option):
            if(not self.continuum_simulating_option):
                self.case_option = 0
            else:
                self.case_option = 3
                
        self.fixed_vel_top = Param.LoadingVelocityTop
        self.fixed_vel_bot = Param.LoadingVelocityBot

        # BOUNDING_BOX
        self.enlargement_factor = Param.BoundingBoxEnlargementFactor
        self.top_corner = Array3()
        self.bottom_corner = Array3()
        self.top_corner[0] = Param.BoundingBoxMaxX
        self.top_corner[0] = Param.BoundingBoxMaxY
        self.top_corner[0] = Param.BoundingBoxMaxZ
        self.bottom_corner[0] = Param.BoundingBoxMinX
        self.bottom_corner[0] = Param.BoundingBoxMinY
        self.bottom_corner[0] = Param.BoundingBoxMinZ

        # MODEL
        self.model_part = model_part
        self.fem_model_part = fem_model_part
        self.contact_model_part = ModelPart("ContactModelPart")  # funcio kratos
        # self.contact_model_part.Nodes       = self.model_part.Nodes; #This is not necessary, elements point at the nodes of the balls_model_part already. It is also problematic when summing modelparts!
        self.domain_size = Param.Dimension


        # GLOBAL PHISICAL ASPECTS
        self.gravity = Vector(3)
        self.gravity[0] = Param.GravityX
        self.gravity[1] = Param.GravityY
        self.gravity[2] = Param.GravityZ

        if (Param.NormalForceCalculationType == "Linear"):
            self.force_calculation_type_id = 0
        elif (Param.NormalForceCalculationType == "Hertz"):
            self.force_calculation_type_id = 1
        elif (Param.NormalForceCalculationType == "PlasticityAndDamage1D"):
            self.force_calculation_type_id = 2
        elif (Param.NormalForceCalculationType == "NonLinearDonze"):
            self.force_calculation_type_id = 3

        if (self.force_calculation_type_id == 2):
            self.C1 = Param.C1
            self.C2 = Param.C2
            self.C3 = Param.C3
            self.N1 = Param.N1
            self.N2 = Param.N2
            self.N3 = Param.N3
            self.plastic_young_modulus_ratio = Param.PlasticYoungModulus
            self.plastic_yield_stress = Param.PlasticYieldStress
            self.damage_deformation_factor = Param.DamageDeformationFactor

        if (self.force_calculation_type_id == 3):
            self.donze_g1 = Param.G1
            self.donze_g2 = Param.G2
            self.donze_g3 = Param.G3
            self.donze_max_def = Param.MaxDef

        if (Param.NormalDampingType == "ViscDamp"):

            if (Param.TangentialDampingType == "ViscDamp"):
                self.damp_id = 11

            else:
                self.damp_id = 10
        else:

            if (Param.TangentialDampingType == "ViscDamp"):
                self.damp_id = 1

            else:
                self.damp_id = 0

        if (Var_Translator(Param.RollingFrictionOption)):
            self.rolling_friction_option = 1

        if (Param.FailureCriterionType == "Mohr-Coulomb"):
            self.failure_criterion_option = 1

        elif (Param.FailureCriterionType == "Uncoupled"):
            self.failure_criterion_option = 2

        self.tau_zero = Param.TauZero
        self.sigma_min = Param.SigmaMin
        self.internal_fricc = Param.InternalFriction

        # PRINTING VARIABLES
        self.print_export_id = Var_Translator(Param.PostExportId)
        self.print_export_skin_sphere = Var_Translator(Param.PostExportSkinSphere)
        self.print_group_id = Var_Translator(Param.PostGroupId)

        self.dummy_switch = 0

        # TIME RELATED PARAMETERS

        self.delta_time = Param.MaxTimeStep
        self.max_delta_time = Param.MaxTimeStep
        self.final_time = Param.FinalTime

        # RESOLUTION METHODS AND PARAMETERS
        self.n_step_search = int(Param.TimeStepsPerSearchStep)
        self.safety_factor = Param.DeltaTimeSafetyFactor  # For critical time step

        # CREATOR-DESTRUCTOR
        self.creator_destructor = creator_destructor

        b_box_low = Array3()
        b_box_high = Array3()
        b_box_low[0] = Param.BoundingBoxMinX
        b_box_low[1] = Param.BoundingBoxMinY
        b_box_low[2] = Param.BoundingBoxMinZ
        b_box_high[0] = Param.BoundingBoxMaxX
        b_box_high[1] = Param.BoundingBoxMaxY
        b_box_high[2] = Param.BoundingBoxMaxZ

        self.creator_destructor.SetLowNode(b_box_low)
        self.creator_destructor.SetHighNode(b_box_high)

        if (self.automatic_bounding_box_option):
            self.creator_destructor.CalculateSurroundingBoundingBox(self.model_part, self.enlargement_factor)

        # STRATEGIES

        self.search_strategy = OMP_DEMSearch()

        if (Param.IntegrationScheme == 'forward_euler'):
            self.time_scheme = ForwardEulerScheme()
        elif (Param.IntegrationScheme == 'mid_point_rule'):
            self.time_scheme = MidPointScheme()
        elif (Param.IntegrationScheme == 'const_average_acc'):
            self.time_scheme = ConstAverageAccelerationScheme()
        else:
            print('scheme not defined')

    #

    def Initialize(self):

        # Setting ProcessInfo variables

        # SIMULATION FLAGS
        self.model_part.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, self.virtual_mass_option)
        self.model_part.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, self.critical_time_option)
        self.model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_option)
        self.model_part.ProcessInfo.SetValue(TRIHEDRON_OPTION, self.trihedron_option)
        self.model_part.ProcessInfo.SetValue(ROTATION_OPTION, self.rotation_option)
        self.model_part.ProcessInfo.SetValue(BOUNDING_BOX_OPTION, self.bounding_box_option)
        self.model_part.ProcessInfo.SetValue(DEMPACK_OPTION, self.dempack_option)
        self.model_part.ProcessInfo.SetValue(ACTIVATE_SEARCH, self.activate_search)
        self.model_part.ProcessInfo.SetValue(FIX_VELOCITIES_FLAG, self.fix_velocities_flag)
        self.model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0);
        self.model_part.ProcessInfo.SetValue(TOTAL_CONTACTS, 0);
        self.model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_option);
       
        # TOLERANCES
        self.model_part.ProcessInfo.SetValue(DISTANCE_TOLERANCE, 0);

# BOUNDARY
#        self.model_part.ProcessInfo.SetValue(LIMIT_SURFACE_OPTION, self.limit_surface_option)
#        if (self.limit_surface_option > 0):
#          self.model_part.ProcessInfo.SetValue(SURFACE_NORMAL_DIR_1, self.surface_normal_dir_1)
#          self.model_part.ProcessInfo.SetValue(SURFACE_POINT_COOR_1, self.surface_point_coor_1)
#          self.model_part.ProcessInfo.SetValue(SURFACE_FRICTION_1, self.surface_friction_angle_1)
#        if (self.limit_surface_option > 1):
#          self.model_part.ProcessInfo.SetValue(SURFACE_NORMAL_DIR_2, self.surface_normal_dir_2)
#          self.model_part.ProcessInfo.SetValue(SURFACE_POINT_COOR_2, self.surface_point_coor_2)
#          self.model_part.ProcessInfo.SetValue(SURFACE_FRICTION_2, self.surface_friction_angle_2)
#        if (self.limit_surface_option > 2):
#          self.model_part.ProcessInfo.SetValue(SURFACE_NORMAL_DIR_3, self.surface_normal_dir_3)
#          self.model_part.ProcessInfo.SetValue(SURFACE_POINT_COOR_3, self.surface_point_coor_3)
#          self.model_part.ProcessInfo.SetValue(SURFACE_FRICTION_3, self.surface_friction_angle_3)
#        if (self.limit_surface_option > 3):
#          self.model_part.ProcessInfo.SetValue(SURFACE_NORMAL_DIR_4, self.surface_normal_dir_4)
#          self.model_part.ProcessInfo.SetValue(SURFACE_POINT_COOR_4, self.surface_point_coor_4)
#          self.model_part.ProcessInfo.SetValue(SURFACE_FRICTION_4, self.surface_friction_angle_4)
#        if (self.limit_surface_option > 4):
#          self.model_part.ProcessInfo.SetValue(SURFACE_NORMAL_DIR_5, self.surface_normal_dir_5)
#          self.model_part.ProcessInfo.SetValue(SURFACE_POINT_COOR_5, self.surface_point_coor_5)
#          self.model_part.ProcessInfo.SetValue(SURFACE_FRICTION_5, self.surface_friction_angle_5)
#
#        self.model_part.ProcessInfo.SetValue(LIMIT_CYLINDER_OPTION, self.limit_cylinder_option)
#        if (self.limit_cylinder_option > 0):
#          self.model_part.ProcessInfo.SetValue(CYLINDER_AXIS_DIR_1, self.cylinder_axis_dir_1)
#          self.model_part.ProcessInfo.SetValue(INITIAL_BASE_CYLINDER_CENTRE_1, self.cylinder_initial_base_centre_1)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_RADIUS_1, self.cylinder_radius_1)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_VELOCITY_1, self.cylinder_velocity_1)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_ANGULAR_VELOCITY_1, self.cylinder_angular_velocity_1)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_FRICTION_1, self.cylinder_friction_angle_1)
#        if (self.limit_cylinder_option > 1):
#          self.model_part.ProcessInfo.SetValue(CYLINDER_AXIS_DIR_2, self.cylinder_axis_dir_2)
#          self.model_part.ProcessInfo.SetValue(INITIAL_BASE_CYLINDER_CENTRE_2, self.cylinder_initial_base_centre_2)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_RADIUS_2, self.cylinder_radius_2)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_VELOCITY_2, self.cylinder_velocity_2)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_ANGULAR_VELOCITY_2, self.cylinder_angular_velocity_2)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_FRICTION_2, self.cylinder_friction_angle_2)
#        if (self.limit_cylinder_option > 2):
#          self.model_part.ProcessInfo.SetValue(CYLINDER_AXIS_DIR_3, self.cylinder_axis_dir_3)
#          self.model_part.ProcessInfo.SetValue(INITIAL_BASE_CYLINDER_CENTRE_3, self.cylinder_initial_base_centre_3)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_RADIUS_3, self.cylinder_radius_3)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_VELOCITY_3, self.cylinder_velocity_3)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_ANGULAR_VELOCITY_3, self.cylinder_angular_velocity_3)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_FRICTION_3, self.cylinder_friction_angle_3)
#        if (self.limit_cylinder_option > 3):
#          self.model_part.ProcessInfo.SetValue(CYLINDER_AXIS_DIR_4, self.cylinder_axis_dir_4)
#          self.model_part.ProcessInfo.SetValue(INITIAL_BASE_CYLINDER_CENTRE_4, self.cylinder_initial_base_centre_4)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_RADIUS_4, self.cylinder_radius_4)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_VELOCITY_4, self.cylinder_velocity_4)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_ANGULAR_VELOCITY_4, self.cylinder_angular_velocity_4)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_FRICTION_4, self.cylinder_friction_angle_4)
#        if (self.limit_cylinder_option > 4):
#          self.model_part.ProcessInfo.SetValue(CYLINDER_AXIS_DIR_5, self.cylinder_axis_dir_5)
#          self.model_part.ProcessInfo.SetValue(INITIAL_BASE_CYLINDER_CENTRE_5, self.cylinder_initial_base_centre_5)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_RADIUS_5, self.cylinder_radius_5)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_VELOCITY_5, self.cylinder_velocity_5)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_ANGULAR_VELOCITY_5, self.cylinder_angular_velocity_5)
#          self.model_part.ProcessInfo.SetValue(CYLINDER_FRICTION_5, self.cylinder_friction_angle_5)
#
        # GLOBAL PHISICAL ASPECTS
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)

        # GLOBAL MATERIAL PROPERTIES

        self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)

        # PRINTING VARIABLES
        self.model_part.ProcessInfo.SetValue(PRINT_GROUP_ID, self.print_group_id)
        self.model_part.ProcessInfo.SetValue(PRINT_EXPORT_ID, self.print_export_id)
        self.model_part.ProcessInfo.SetValue(PRINT_SKIN_SPHERE, self.print_export_skin_sphere)
        self.model_part.ProcessInfo.SetValue(FORCE_CALCULATION_TYPE, self.force_calculation_type_id)
        self.model_part.ProcessInfo.SetValue(DAMP_TYPE, self.damp_id)
        self.model_part.ProcessInfo.SetValue(ROLLING_FRICTION_OPTION, self.rolling_friction_option)

        # TIME RELATED PARAMETERS
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)
        self.model_part.ProcessInfo.SetValue(FINAL_SIMULATION_TIME, self.final_time)

        # CONTINUUM

        self.model_part.ProcessInfo.SetValue(SEARCH_TOLERANCE, self.search_tolerance)
        self.model_part.ProcessInfo.SetValue(AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION, self.amplified_continuum_search_radius_extension)

        self.model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, self.contact_mesh_option)

        self.model_part.ProcessInfo.SetValue(FAILURE_CRITERION_OPTION, self.failure_criterion_option)
        self.model_part.ProcessInfo.SetValue(CONTACT_SIGMA_MIN, self.sigma_min)
        self.model_part.ProcessInfo.SetValue(CONTACT_TAU_ZERO, self.tau_zero)
        self.model_part.ProcessInfo.SetValue(CONTACT_INTERNAL_FRICC, self.internal_fricc)

        if(self.dempack_option):
            self.model_part.ProcessInfo.SetValue(DEMPACK_DAMPING, self.dempack_damping)
            self.model_part.ProcessInfo.SetValue(DEMPACK_GLOBAL_DAMPING, self.dempack_global_damping)
            self.model_part.ProcessInfo.SetValue(AREA_FACTOR, self.area_factor)

        if (self.force_calculation_type_id == 2):
            self.model_part.ProcessInfo.SetValue(SLOPE_FRACTION_N1, self.N1)
            self.model_part.ProcessInfo.SetValue(SLOPE_FRACTION_N2, self.N2)
            self.model_part.ProcessInfo.SetValue(SLOPE_FRACTION_N3, self.N3)
            self.model_part.ProcessInfo.SetValue(SLOPE_LIMIT_COEFF_C1, self.C1)
            self.model_part.ProcessInfo.SetValue(SLOPE_LIMIT_COEFF_C2, self.C2)
            self.model_part.ProcessInfo.SetValue(SLOPE_LIMIT_COEFF_C3, self.C3)
            self.model_part.ProcessInfo.SetValue(YOUNG_MODULUS_PLASTIC, self.plastic_young_modulus_ratio)
            self.model_part.ProcessInfo.SetValue(PLASTIC_YIELD_STRESS, self.plastic_yield_stress)
            self.model_part.ProcessInfo.SetValue(DAMAGE_FACTOR, self.damage_deformation_factor)

        if (self.force_calculation_type_id == 3):
            self.model_part.ProcessInfo.SetValue(DONZE_G1, self.donze_g1)
            self.model_part.ProcessInfo.SetValue(DONZE_G2, self.donze_g2)
            self.model_part.ProcessInfo.SetValue(DONZE_G3, self.donze_g3)
            self.model_part.ProcessInfo.SetValue(DONZE_MAX_DEF, self.donze_max_def)

        if ( (self.concrete_test_option == "TRIAXIAL") or (self.concrete_test_option == "HYDROSTATIC")):
            self.model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 1)
            
        self.model_part.ProcessInfo.SetValue(FIXED_VEL_TOP, self.fixed_vel_top)
        self.model_part.ProcessInfo.SetValue(FIXED_VEL_BOT, self.fixed_vel_bot)
            

        # OTHERS

        self.model_part.ProcessInfo.SetValue(DUMMY_SWITCH, self.dummy_switch)

        # RESOLUTION METHODS AND PARAMETERS
        # Creating the solution strategy

        self.solver = ContinuumExplicitSolverStrategy(self.model_part, self.fem_model_part, self.contact_model_part, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                      self.MoveMeshFlag, self.delta_option, self.search_tolerance, self.coordination_number, self.creator_destructor, self.time_scheme, self.search_strategy)

        self.solver.Initialize()  # Calls the solver Initialized function (initializes all elements and performs other necessary tasks before iterating)

    #

    def Initial_Critical_Time(self):
        (self.solver).InitialTimeStepCalculation()

    #

    def Solve(self):
        (self.solver).Solve()
