# Model Type

ContinuumOption= "*GenData(Continuum)"
RotationOption= "*GenData(Rotation)"
HomogeneousMaterial= "*GenData(Homogeneous_Material)"

# General Settings

*format "%10.5e"
final_time = *GenData(Calculation_Time)
*format "%10.5e"
gravity_x = *GenData(Gravity_x)
*format "%10.5e"
gravity_y = *GenData(Gravity_y)
*format "%10.5e"
gravity_z = *GenData(Gravity_z)
BoundingBoxOption= "*GenData(Bounding_Box_option)"
*format "%10.5e"
bounding_box_enlargement_factor = *GenData(Bounding_Box_Enlargement_Factor)
domain_size = *GenData(Dimensions)
OutputFileType="*GenData(Output_file_type)"
Multifile="*GenData(Multifile)"
PrintNeighbourLists= "*GenData(Print_Neighbour_Lists)"
ModelDataInfo = "*GenData(Model_Data_Info)"

# Special features

VirtualMassOption = "*GenData(Virtual_Mass_Option_Id)"
*format "%10.5e"
VirtualMassCoefficient = *GenData(Virtual_Mass_Coefficient)
*format "%10.5e"
MagicFactor = *GenData(Magic_Factor)
DeltaOption = "*GenData(Delta_Option)"
*format "%10.5e"
search_radius_extension = *GenData(Search_Radius_Extension)
FixVelocities = "*GenData(Fix_Velocities_At_Predetermined_Time)"
*format "%10.5e"
TimePercentageFixVelocities = *GenData(Time_Step_Constrain_DOFs_Percentage)
TrihedronOption= "*GenData(Trihedron_Option)"
LimitSurfaceOption= "*GenData(Limit_Surface_Option)"
*format "%10.5e"
surface_normal_dir_x = *GenData(Surface_Normal_Dir_X)
*format "%10.5e"
surface_normal_dir_y = *GenData(Surface_Normal_Dir_Y)
*format "%10.5e"
surface_normal_dir_z = *GenData(Surface_Normal_Dir_Z)
*format "%10.5e"
surface_point_coor_x = *GenData(Surface_Point_Coor_X)
*format "%10.5e"
surface_point_coor_y = *GenData(Surface_Point_Coor_Y)
*format "%10.5e"
surface_point_coor_z = *GenData(Surface_Point_Coor_Z)
surface_friction_angle = *GenData(Surface_Friction_Angle)

# Time Discretization Settings

Integration_Scheme = "*GenData(Integration_Scheme)"
*format "%10.5e"
search_step = *GenData(Search_Step)
CriticalTimeOption = "*GenData(Critical_Time_Option)"
*format "%10.5e"
dt_safety_factor = *GenData(Dt_Safety_Factor)
*format "%10.5e"
max_time_step = *GenData(Max_Time_Step)
*format "%10.5e"
output_dt = *GenData(Output_Dt)
*format "%10.5e"
control_time = *GenData(Control_Time)

# Material Model

NormalForceCalculation = "*GenData(Normal_Force_Calculation)"
NormalDampId = "*GenData(Normal_Contact_Damp)"
TangentialDampId = "*GenData(Tangential_Contact_Damp)"
FailureCriterionOption = "*GenData(Failure_Criterion_Option)"
*format "%10.5e"
TauZero = *GenData(Tau_Zero)
*format "%10.5e"
SigmaMax = *GenData(Sigma_Max)
*format "%10.5e"
SigmaMin = *GenData(Sigma_Min)
*format "%10.5e"
InternalFricc = *GenData(Internal_Fricc)
NonLinearOption = "*GenData(Non_Linear_Option)"
*format "%10.5e"
C1 = *GenData(C_1)
*format "%10.5e"
N1 = *GenData(N_1)
*format "%10.5e"
C2 = *GenData(C_2)
*format "%10.5e"
N2 = *GenData(N_2)
RotationalSpringOption= "*GenData(Rotation_Spring)"
RotaDampId = "*GenData(Rota_Damp_Id)"

# Global Material Parameters

*format "%10.5e"
GeneralDensity= "*GenData(General_Density)"
*format "%10.5e"
GeneralYoungModulus= "*GenData(General_Young_Modulus)"
*format "%10.5e"
GeneralPoissonRatio= "*GenData(General_Poisson_Ratio)"
*format "%10.5e"
GeneralCohesion= "*GenData(General_Cohesion)"
*format "%10.5e"
GeneralRollingFriction= "*GenData(General_Rolling_Friction)"
*format "%10.5e"
GeneralTension= "*GenData(General_Tension)"
*format "%10.5e"
GeneralRotaDampRatio= "*GenData(General_Rota_Damp_Ratio)"
*format "%10.5e"
GeneralStaticFrictionCoef= "*GenData(General_Static_Friction_Coef)"
*format "%10.5e"
GeneralDynamicFrictionCoef= "*GenData(General_Dynamic_Friction_Coef)"
*format "%10.5e"
GeneralRestitutionCoef= "*GenData(General_Restitution_Coef)"
*format "%10.5e"
GeneralColour= "*GenData(General_Colour)"
*format "%10.5e"
GlobalVariablesOption= "*GenData(Globally_Specified_Variables)"
*format "%10.5e"
global_kn = *GenData(Global_KN)
*format "%10.5e"
global_kt = *GenData(Global_KT)
*format "%10.5e"
global_kr = *GenData(Global_KR)
*format "%10.5e"
global_rn = *GenData(Global_RN)
*format "%10.5e"
global_rt = *GenData(Global_RT)
*format "%10.5e"
global_rr = *GenData(Global_RR)
*format "%10.5e"
global_fri_ang = *GenData(Global_FRI_ANG)

# Continuum Options

StressStrainOperations= "*GenData(Stress_Strain_Operations)"
ContactMeshOption = "*GenData(Contact_Mesh_Option)"
ConcreteTestOption = "*GenData(Concrete_Test_Option)"
RealTimeGraph = "*GenData(Real_Time_Graph)"
TriaxialOption = "*GenData(Triaxial_Option)"
*format "%10.5e"
ConfinementPressure = *GenData(Confinement_Pressure)
*format "%10.5e"
InitialTime = *GenData(Initial_Time)
*format "%10.5e"
IncreasingTemporaily = *GenData(Increasing_Temporaily)

#POSTPROCES

print_velocity                   = "*GenData(VELOCITY)"
print_displacement               = "*GenData(DISPLACEMENT)"
print_radial_displacement        = "*GenData(RADIAL_DISPLACEMENT)"
print_rhs                        = "*GenData(RHS)"
print_total_forces               = "*GenData(TOTAL_FORCES)"
print_damp_forces                = "*GenData(DAMP_FORCES)"
print_applied_forces             = "*GenData(APPLIED_FORCES)"
print_radius                     = "*GenData(RADIUS)"
print_particle_cohesion          = "*GenData(PARTICLE_COHESION)"
print_particle_tension           = "*GenData(PARTICLE_TENSION)"
print_group_id                   = "*GenData(GROUP_ID)"
print_export_id                  = "*GenData(EXPORT_ID)"
print_export_particle_failure_id = "*GenData(EXPORT_PARTICLE_FAILURE_ID)"
print_export_skin_sphere         = "*GenData(EXPORT_SKIN_SPHERE)"
print_local_contact_force_low    = "*GenData(LOCAL_CONTACT_FORCE_LOW)"
print_local_contact_force_high   = "*GenData(LOCAL_CONTACT_FORCE_HIGH)"
print_failure_criterion_state    = "*GenData(FAILURE_CRITERION_STATE)"
print_contact_failure            = "*GenData(CONTACT_FAILURE)"
print_contact_tau                = "*GenData(CONTACT_TAU)"
print_contact_sigma              = "*GenData(CONTACT_SIGMA)"
print_angular_velocity           = "*GenData(ANGULAR_VELOCITY)"
print_particle_moment            = "*GenData(PARTICLE_MOMENT)"
print_euler_angles               = "*GenData(EULER_ANGLES)"
print_representative_volume      = "*GenData(REPRESENTATIVE_VOLUME)"
print_mean_contact_area          = "*GenData(MEAN_CONTACT_AREA)"
print_stress_tensor              = "*GenData(STRESS_TENSOR)"

#FROM CND:

*Set cond SET_SKIN *elems
*if(CondNumEntities(int))
predefined_skin_option = "ON"
*else
predefined_skin_option = "OFF"
*endif

mass_elements=*tcl(DEM::Get_Mass_Elements)

# Declare Python Variables
