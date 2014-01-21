from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *


def Var_Translator(variable):

    if (variable == "OFF" or variable == "0" or variable == 0):
        variable = 0
    else:
        variable = 1

    return variable


def AddVariables(model_part, Param):

    # KINEMATIC
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(RADIAL_DISPLACEMENT)
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
        model_part.AddNodalSolutionStepVariable(PARTICLE_INERTIA)
        model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA)
        model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_DAMP_RATIO)
        model_part.AddNodalSolutionStepVariable(ROLLING_FRICTION)

    # OTHER PROPERTIES
    model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD
    model_part.AddNodalSolutionStepVariable(PARTICLE_CONTINUUM)  # Continuum group
#    model_part.AddNodalSolutionStepVariable(REPRESENTATIVE_VOLUME)
    model_part.AddNodalSolutionStepVariable(MAX_INDENTATION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_SPHERICITY)  # MA: this is added temporarily until inlet becomes a process

    # LOCAL AXIS
    model_part.AddNodalSolutionStepVariable(EULER_ANGLES)

    # BOUNDARY SURFACE

#   if (Param.LimitSurfaceOption > 0):
#       model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES_1)
#   if (Param.LimitSurfaceOption > 1):
#       model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES_2)
#   if (Param.LimitSurfaceOption > 2):
#       model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES_3)
#   if (Param.LimitSurfaceOption > 3):
#       model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES_4)
#   if (Param.LimitSurfaceOption > 4):
#       model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES_5)
#
#   if (Param.LimitCylinderOption > 0):
#       model_part.AddNodalSolutionStepVariable(PARTICLE_CYLINDER_CONTACT_FORCES_1)
#   if (Param.LimitCylinderOption > 1):
#       model_part.AddNodalSolutionStepVariable(PARTICLE_CYLINDER_CONTACT_FORCES_2)
#   if (Param.LimitCylinderOption > 2):
#       model_part.AddNodalSolutionStepVariable(PARTICLE_CYLINDER_CONTACT_FORCES_3)
#   if (Param.LimitCylinderOption > 3):
#       model_part.AddNodalSolutionStepVariable(PARTICLE_CYLINDER_CONTACT_FORCES_4)
#   if (Param.LimitCylinderOption > 4):
#       model_part.AddNodalSolutionStepVariable(PARTICLE_CYLINDER_CONTACT_FORCES_5)

    # OPTIMIZATION
    model_part.AddNodalSolutionStepVariable(VELOCITY_X_DOF_POS)
    model_part.AddNodalSolutionStepVariable(VELOCITY_Y_DOF_POS)
    model_part.AddNodalSolutionStepVariable(VELOCITY_Z_DOF_POS)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY_X_DOF_POS)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY_Y_DOF_POS)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY_Z_DOF_POS)
    model_part.AddNodalSolutionStepVariable(OLD_COORDINATES)

    # FLAGS
    model_part.AddNodalSolutionStepVariable(GROUP_ID)            # Differencied groups for plotting, etc..
#    model_part.AddNodalSolutionStepVariable(ERASE_FLAG)

    # ONLY VISUALIZATION
    model_part.AddNodalSolutionStepVariable(EXPORT_ID)

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
        node.AddDof(ANGULAR_VELOCITY_X, REACTION_X)
        node.AddDof(ANGULAR_VELOCITY_Y, REACTION_Y)
        node.AddDof(ANGULAR_VELOCITY_Z, REACTION_Z)

    print("DOFs for the DEM solution added correctly")


class ExplicitStrategy:

    def __init__(self, model_part, fem_model_part, creator_destructor, Param):

        # Initialization of member variables

        # SIMULATION FLAGS
        self.virtual_mass_option = Var_Translator(Param.VirtualMassOption)
        self.critical_time_option = Var_Translator(Param.AutoReductionOfTimeStepOption)
        self.trihedron_option = Var_Translator(Param.PostEulerAngles)
        self.rotation_option = Var_Translator(Param.RotationOption)
        self.bounding_box_option = Var_Translator(Param.BoundingBoxOption)
        self.fix_velocities = Var_Translator(Param.FixVelocitiesOption)
#        self.limit_surface_option           = Param.LimitSurfaceOption
#        self.limit_cylinder_option          = Param.LimitCylinderOption
        self.clean_init_indentation_option = Var_Translator(Param.CleanIndentationsOption)
        self.homogeneous_material_option = Var_Translator(Param.HomogeneousMaterialOption)
        self.global_variables_option = Var_Translator(Param.GlobalVariablesOption)
        self.contact_mesh_option = Var_Translator(Param.ContactMeshOption)
        self.automatic_bounding_box_option = Var_Translator(Param.AutomaticBoundingBoxOption)

        if (Param.DeltaOption == "OFF"):
            self.delta_option = 0

        elif (Param.DeltaOption == "ABSOLUTE"):
            self.delta_option = 1
            self.search_tolerance = Param.SearchTolerance

        elif (Param.DeltaOption == "COORDINATION_NUMBER"):
            self.delta_option = 2
            self.coordination_number = Param.CoordinationNumber

        self.move_mesh_flag = True
        self.deactivate_search = 0
        self.case_option = 3

        self.ComputeMovementOption = Var_Translator(Param.ComputeMovementOption)

        # MODEL
        self.model_part = model_part
        self.fem_model_part = fem_model_part

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

        # RigidFace movement
        if(self.ComputeMovementOption):
            self.RotationSpeed = Param.RotationSpeed
            self.AxialSpeed = Param.AxialSpeed
            self.PROP_ID = Param.PropID
            self.GLOBAL_VEL = Vector(3)
            self.GLOBAL_VEL[0] = Param.GLOBAL_X_VEL
            self.GLOBAL_VEL[1] = Param.GLOBAL_Y_VEL
            self.GLOBAL_VEL[2] = Param.GLOBAL_Z_VEL
            self.ROTA_ORIGIN_COORD = Vector(3)
            self.ROTA_ORIGIN_COORD[0] = Param.ROTA_ORIGIN_COORD_X
            self.ROTA_ORIGIN_COORD[1] = Param.ROTA_ORIGIN_COORD_Y
            self.ROTA_ORIGIN_COORD[2] = Param.ROTA_ORIGIN_COORD_Z
            self.ROTA_AXIAL_NORMAL = Vector(3)
            self.ROTA_AXIAL_NORMAL[0] = Param.ROTA_AXIAL_NORMAL_X
            self.ROTA_AXIAL_NORMAL[1] = Param.ROTA_AXIAL_NORMAL_Y
            self.ROTA_AXIAL_NORMAL[2] = Param.ROTA_AXIAL_NORMAL_Z
            self.BEGIN_TIME = Param.BEGIN_TIME
            self.END_TIME = Param.END_TIME

        # BOUNDARY
#        if (Param.LimitSurfaceOption > 0):
#          self.surface_normal_dir_1           = Vector(3)
#          self.surface_normal_dir_1[0]        = Param.SurfaceNormalDirX1
#          self.surface_normal_dir_1[1]        = Param.SurfaceNormalDirY1
#          self.surface_normal_dir_1[2]        = Param.SurfaceNormalDirZ1
#          self.surface_point_coor_1           = Vector(3)
#          self.surface_point_coor_1[0]        = Param.SurfacePointCoorX1
#          self.surface_point_coor_1[1]        = Param.SurfacePointCoorY1
#          self.surface_point_coor_1[2]        = Param.SurfacePointCoorZ1
#          self.surface_friction_angle_1       = Param.SurfaceFrictionAngle1
#        if (Param.LimitSurfaceOption > 1):
#          self.surface_normal_dir_2           = Vector(3)
#          self.surface_normal_dir_2[0]        = Param.SurfaceNormalDirX2
#          self.surface_normal_dir_2[1]        = Param.SurfaceNormalDirY2
#          self.surface_normal_dir_2[2]        = Param.SurfaceNormalDirZ2
#          self.surface_point_coor_2           = Vector(3)
#          self.surface_point_coor_2[0]        = Param.SurfacePointCoorX2
#          self.surface_point_coor_2[1]        = Param.SurfacePointCoorY2
#          self.surface_point_coor_2[2]        = Param.SurfacePointCoorZ2
#          self.surface_friction_angle_2       = Param.SurfaceFrictionAngle2
#        if (Param.LimitSurfaceOption > 2):
#          self.surface_normal_dir_3           = Vector(3)
#          self.surface_normal_dir_3[0]        = Param.SurfaceNormalDirX3
#          self.surface_normal_dir_3[1]        = Param.SurfaceNormalDirY3
#          self.surface_normal_dir_3[2]        = Param.SurfaceNormalDirZ3
#          self.surface_point_coor_3           = Vector(3)
#          self.surface_point_coor_3[0]        = Param.SurfacePointCoorX3
#          self.surface_point_coor_3[1]        = Param.SurfacePointCoorY3
#          self.surface_point_coor_3[2]        = Param.SurfacePointCoorZ3
#          self.surface_friction_angle_3       = Param.SurfaceFrictionAngle3
#        if (Param.LimitSurfaceOption > 3):
#          self.surface_normal_dir_4           = Vector(3)
#          self.surface_normal_dir_4[0]        = Param.SurfaceNormalDirX4
#          self.surface_normal_dir_4[1]        = Param.SurfaceNormalDirY4
#          self.surface_normal_dir_4[2]        = Param.SurfaceNormalDirZ4
#          self.surface_point_coor_4           = Vector(3)
#          self.surface_point_coor_4[0]        = Param.SurfacePointCoorX4
#          self.surface_point_coor_4[1]        = Param.SurfacePointCoorY4
#          self.surface_point_coor_4[2]        = Param.SurfacePointCoorZ4
#          self.surface_friction_angle_4       = Param.SurfaceFrictionAngle4
#        if (Param.LimitSurfaceOption > 4):
#          self.surface_normal_dir_5           = Vector(3)
#          self.surface_normal_dir_5[0]        = Param.SurfaceNormalDirX5
#          self.surface_normal_dir_5[1]        = Param.SurfaceNormalDirY5
#          self.surface_normal_dir_5[2]        = Param.SurfaceNormalDirZ5
#          self.surface_point_coor_5           = Vector(3)
#          self.surface_point_coor_5[0]        = Param.SurfacePointCoorX5
#          self.surface_point_coor_5[1]        = Param.SurfacePointCoorY5
#          self.surface_point_coor_5[2]        = Param.SurfacePointCoorZ5
#          self.surface_friction_angle_5       = Param.SurfaceFrictionAngle5
#
#        if (Param.LimitCylinderOption > 0):
#          self.cylinder_axis_dir_1              = Vector(3)
#          self.cylinder_axis_dir_1[0]           = Param.CylinderAxisX1
#          self.cylinder_axis_dir_1[1]           = Param.CylinderAxisY1
#          self.cylinder_axis_dir_1[2]           = Param.CylinderAxisZ1
#          self.cylinder_initial_base_centre_1   = Vector(3)
#          self.cylinder_initial_base_centre_1[0]= Param.CylinderInitialBaseCentreX1
#          self.cylinder_initial_base_centre_1[1]= Param.CylinderInitialBaseCentreY1
#          self.cylinder_initial_base_centre_1[2]= Param.CylinderInitialBaseCentreZ1
#          self.cylinder_radius_1                = Param.CylinderRadius1
#          self.cylinder_velocity_1              = Param.CylinderVelocity1
#          self.cylinder_angular_velocity_1      = Param.CylinderAngularVelocity1
#          self.cylinder_friction_angle_1        = Param.CylinderFrictionAngle1
#        if (Param.LimitCylinderOption > 1):
#          self.cylinder_axis_dir_2              = Vector(3)
#          self.cylinder_axis_dir_2[0]           = Param.CylinderAxisX2
#          self.cylinder_axis_dir_2[1]           = Param.CylinderAxisY2
#          self.cylinder_axis_dir_2[2]           = Param.CylinderAxisZ2
#          self.cylinder_initial_base_centre_2   = Vector(3)
#          self.cylinder_initial_base_centre_2[0]= Param.CylinderInitialBaseCentreX2
#          self.cylinder_initial_base_centre_2[1]= Param.CylinderInitialBaseCentreY2
#          self.cylinder_initial_base_centre_2[2]= Param.CylinderInitialBaseCentreZ2
#          self.cylinder_radius_2                = Param.CylinderRadius2
#          self.cylinder_velocity_2              = Param.CylinderVelocity2
#          self.cylinder_angular_velocity_2      = Param.CylinderAngularVelocity2
#          self.cylinder_friction_angle_2        = Param.CylinderFrictionAngle2
#        if (Param.LimitCylinderOption > 2):
#          self.cylinder_axis_dir_3              = Vector(3)
#          self.cylinder_axis_dir_3[0]           = Param.CylinderAxisX3
#          self.cylinder_axis_dir_3[1]           = Param.CylinderAxisY3
#          self.cylinder_axis_dir_3[2]           = Param.CylinderAxisZ3
#          self.cylinder_initial_base_centre_3   = Vector(3)
#          self.cylinder_initial_base_centre_3[0]= Param.CylinderInitialBaseCentreX3
#          self.cylinder_initial_base_centre_3[1]= Param.CylinderInitialBaseCentreY3
#          self.cylinder_initial_base_centre_3[2]= Param.CylinderInitialBaseCentreZ3
#          self.cylinder_radius_3                = Param.CylinderRadius3
#          self.cylinder_velocity_3              = Param.CylinderVelocity3
#          self.cylinder_angular_velocity_3      = Param.CylinderAngularVelocity3
#          self.cylinder_friction_angle_3        = Param.CylinderFrictionAngle3
#        if (Param.LimitCylinderOption > 3):
#          self.cylinder_axis_dir_4              = Vector(3)
#          self.cylinder_axis_dir_4[0]           = Param.CylinderAxisX4
#          self.cylinder_axis_dir_4[1]           = Param.CylinderAxisY4
#          self.cylinder_axis_dir_4[2]           = Param.CylinderAxisZ4
#          self.cylinder_initial_base_centre_4   = Vector(3)
#          self.cylinder_initial_base_centre_4[0]= Param.CylinderInitialBaseCentreX4
#          self.cylinder_initial_base_centre_4[1]= Param.CylinderInitialBaseCentreY4
#          self.cylinder_initial_base_centre_4[2]= Param.CylinderInitialBaseCentreZ4
#          self.cylinder_radius_4                = Param.CylinderRadius4
#          self.cylinder_velocity_4              = Param.CylinderVelocity4
#          self.cylinder_angular_velocity_4      = Param.CylinderAngularVelocity4
#          self.cylinder_friction_angle_4        = Param.CylinderFrictionAngle4
#        if (Param.LimitCylinderOption > 4):
#          self.cylinder_axis_dir_5              = Vector(3)
#          self.cylinder_axis_dir_5[0]           = Param.CylinderAxisX5
#          self.cylinder_axis_dir_5[1]           = Param.CylinderAxisY5
#          self.cylinder_axis_dir_5[2]           = Param.CylinderAxisZ5
#          self.cylinder_initial_base_centre_5   = Vector(3)
#          self.cylinder_initial_base_centre_5[0]= Param.CylinderInitialBaseCentreX5
#          self.cylinder_initial_base_centre_5[1]= Param.CylinderInitialBaseCentreY5
#          self.cylinder_initial_base_centre_5[2]= Param.CylinderInitialBaseCentreZ5
#          self.cylinder_radius_5                = Param.CylinderRadius5
#          self.cylinder_velocity_5              = Param.CylinderVelocity5
#          self.cylinder_angular_velocity_5      = Param.CylinderAngularVelocity5
#          self.cylinder_friction_angle_5        = Param.CylinderFrictionAngle5

        # GLOBAL PHYSICAL ASPECTS
        self.gravity = Vector(3)
        self.gravity[0] = Param.GravityX
        self.gravity[1] = Param.GravityY
        self.gravity[2] = Param.GravityZ

        # GLOBAL MATERIAL PROPERTIES
        self.nodal_mass_coeff = Param.VirtualMassCoefficient
        self.magic_factor = Param.MagicFactor

        if (self.global_variables_option):
            self.global_kn = Param.GlobalKn
            self.global_kt = Param.GlobalKt

        if (Param.NormalForceCalculationType == "Linear"):
            self.force_calculation_type_id = 0

        elif (Param.NormalForceCalculationType == "Hertz"):
            self.force_calculation_type_id = 1

        else:

            raise 'Specified NormalForceCalculationType is not defined'

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

        if (Param.RotaDampingType == "LocalDamp"):
            self.rota_damp_id = 1

        elif (Param.RotaDampingType == "RollingFric"):
            self.rota_damp_id = 2

        else:
            self.rota_damp_id = 0


        self.tau_zero                       = Param.TauZero

        # PRINTING VARIABLES
        self.print_export_id = Var_Translator(Param.PostExportId)
        self.print_group_id = Var_Translator(Param.PostGroupId)
        self.print_radial_displacement = Var_Translator(Param.PostRadialDisplacement)

        # TIME RELATED PARAMETERS
        self.delta_time = Param.MaxTimeStep
        self.max_delta_time = Param.MaxTimeStep
        self.final_time = Param.FinalTime

        # RESOLUTION METHODS AND PARAMETERS

        if (Param.TimeStepsPerSearchStep < 1):

            raise 'Variable TimeStepsPerSearchStep must be an integer, grater or equal to 1. The current input value is ', Param.TimeStepsPerSearchStep

        elif (not isinstance(Param.TimeStepsPerSearchStep, int)):

            print('Variable TimeStepsPerSearchStep is not an integer. Its input value is ', Param.TimeStepsPerSearchStep, 'Rounding up to ', int(Param.TimeStepsPerSearchStep))

            self.n_step_search = int(Param.TimeStepsPerSearchStep)

        else:
            self.n_step_search = int(Param.TimeStepsPerSearchStep)

        if (self.deactivate_search):
            self.n_step_search = sys.maxsize

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

        self.search_tolerance = 0.0
        self.coordination_number = 10.0

        self.search_strategy = OMP_DEMSearch()

        if (Param.IntegrationScheme == 'forward_euler'):
            self.time_scheme = ForwardEulerScheme()

        elif (Param.IntegrationScheme == 'mid_point_rule'):
            self.time_scheme = MidPointScheme()

        elif (Param.IntegrationScheme == 'const_average_acc'):
            self.time_scheme = ConstAverageAccelerationScheme()

        else:

            print('Specified IntegrationScheme is not defined')

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
        self.model_part.ProcessInfo.SetValue(FIX_VELOCITIES_FLAG, self.fix_velocities)
        self.model_part.ProcessInfo.SetValue(GLOBAL_VARIABLES_OPTION, self.global_variables_option)
        self.model_part.ProcessInfo.SetValue(UNIFORM_MATERIAL_OPTION, self.homogeneous_material_option)
        self.model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0);
        self.model_part.ProcessInfo.SetValue(TOTAL_CONTACTS, 0);
        self.model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_option)
        self.model_part.ProcessInfo.SetValue(RIGID_FACE_FLAG, self.ComputeMovementOption)
        self.model_part.ProcessInfo.SetValue(ACTIVATE_SEARCH, 1)  # needed in the basic for the continuum.

        # TOTAL NUMBER OF INITIALIZED ELEMENTS
        self.model_part.ProcessInfo.SetValue(NUM_PARTICLES_INITIALIZED, 0);

        # TOLERANCES
        self.model_part.ProcessInfo.SetValue(DISTANCE_TOLERANCE, 0);

        # Rigid Face setting
        if(self.ComputeMovementOption):
            self.model_part.ProcessInfo.SetValue(RIGID_FACE_ROTA_SPEED, self.RotationSpeed)
            self.model_part.ProcessInfo.SetValue(RIGID_FACE_AXIAL_SPEED, self.AxialSpeed)
            self.model_part.ProcessInfo.SetValue(RIGID_FACE_PROP_ID, self.PROP_ID)
            self.model_part.ProcessInfo.SetValue(RIGID_FACE_ROTA_ORIGIN_COORD, self.ROTA_ORIGIN_COORD)
            self.model_part.ProcessInfo.SetValue(RIGID_FACE_ROTA_AXIAL_DIR, self.ROTA_AXIAL_NORMAL)
            self.model_part.ProcessInfo.SetValue(RIGID_FACE_ROTA_GLOBAL_VELOCITY, self.GLOBAL_VEL)
            self.model_part.ProcessInfo.SetValue(RIGID_FACE_BEGIN_TIME, self.BEGIN_TIME)
            self.model_part.ProcessInfo.SetValue(RIGID_FACE_END_TIME, self.END_TIME)

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

        # GLOBAL PHISICAL ASPECTS
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)
        self.model_part.ProcessInfo.SetValue(DEM_MAGIC_FACTOR, self.magic_factor)

        # GLOBAL MATERIAL PROPERTIES
        self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)

        if (self.global_variables_option):
            self.model_part.ProcessInfo.SetValue(GLOBAL_KN, self.global_kn)
            self.model_part.ProcessInfo.SetValue(GLOBAL_KT, self.global_kt)

        # SEARCH-RELATED
        self.model_part.ProcessInfo.SetValue(SEARCH_TOLERANCE, self.search_tolerance)  # needed in ProcessInfo for MPISearch

        # PRINTING VARIABLES

        self.model_part.ProcessInfo.SetValue(FORCE_CALCULATION_TYPE, self.force_calculation_type_id)
        self.model_part.ProcessInfo.SetValue(DAMP_TYPE, self.damp_id)
        self.model_part.ProcessInfo.SetValue(ROTA_DAMP_TYPE, self.rota_damp_id)
        self.model_part.ProcessInfo.SetValue(PRINT_RADIAL_DISPLACEMENT, self.print_radial_displacement)
        self.model_part.ProcessInfo.SetValue(PRINT_GROUP_ID, self.print_group_id)
        self.model_part.ProcessInfo.SetValue(PRINT_EXPORT_ID, self.print_export_id)

        # TIME RELATED PARAMETERS
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)
        self.model_part.ProcessInfo.SetValue(FINAL_SIMULATION_TIME, self.final_time)

        # RESOLUTION METHODS AND PARAMETERS
        # Creating the solution strategy

        self.solver = ExplicitSolverStrategy(self.model_part, self.fem_model_part, self.max_delta_time, self.n_step_search, self.safety_factor, self.move_mesh_flag,
                                             self.delta_option, self.search_tolerance, self.coordination_number, self.creator_destructor, self.time_scheme, self.search_strategy)

        self.solver.Initialize()  # Calls the solver Initialize function (initializes all elements and performs other necessary tasks before iterating) (C++)

    #

    def Solve(self):
        (self.solver).Solve()

    def Compute_RigidFace_Movement(self):
        (self.solver).Compute_RigidFace_Movement()

    #
