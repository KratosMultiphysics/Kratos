# General Settings

BoundingBoxOption= "*GenData(Bounding_Box_option)"
*format "%10.5e"
bounding_box_enlargement_factor = *GenData(Bounding_Box_Enlargement_Factor)
domain_size = *GenData(Dimensions)
OutputFileType="*GenData(Output_file_type)"
Multifile="*GenData(Multifile)"

# Calculation Settings

Integration_Scheme = "*GenData(Integration_Scheme)"
NormalForceCalculation = "*GenData(Normal_Force_Calculation)"
DampId = "*GenData(Damp_Id)"
*format "%10.5e"
gravity_x = *GenData(Gravity_x)
*format "%10.5e"
gravity_y = *GenData(Gravity_y)
*format "%10.5e"
gravity_z = *GenData(Gravity_z)

#Time Settings

*format "%10.5e"
search_step = *GenData(Search_Step)
*format "%10.5e"
final_time = *GenData(Calculation_Time)
CriticalTimeOption = "*GenData(Critical_Time_Option)"
*format "%10.5e"
dt_safety_factor = *GenData(Dt_Safety_Factor)
*format "%10.5e"
max_time_step = *GenData(Max_Time_Step)
*format "%10.5e"
output_dt = *GenData(Output_Dt)
*format "%10.5e"
control_time = *GenData(Control_Time)

#Special features

VirtualMassOption = "*GenData(Virtual_Mass_Option_Id)"
*format "%10.5e"
VirtualMassCoefficient = *GenData(Virtual_Mass_Coefficient)
DeltaOption = "*GenData(Delta_Option)"
*format "%10.5e"
search_radius_extension = *GenData(Search_Radius_Extension)
PrintNeighbourLists= "*GenData(Print_Neighbour_Lists)"
GlobalVariablesOption= "*GenData(Global_Variables_Option)"
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
ModelDataInfo = "*GenData(Model_Data_Info)"

#Continuum Options

ContinuumOption= "*GenData(Continuum_Option)"
ContactMeshOption = "*GenData(Contact_Mesh_Option)"
FailureCriterionOption = "*GenData(Failure_Criterion_Option)"
*format "%10.5e"
TauZero = *GenData(Tau_Zero)
*format "%10.5e"
SigmaMax = *GenData(Sigma_Max)
*format "%10.5e"
SigmaMin = *GenData(Sigma_Min)
*format "%10.5e"
InternalFricc = *GenData(Internal_Fricc)

#Concrete Test Options

ConcreteTestOption = "*GenData(Concrete_Test_Option)"
RealTimeGraph = "*GenData(Real_Time_Graph)"
TriaxialOption = "*GenData(Triaxial_Option)"
*format "%10.5e"
ConfinementPressure = *GenData(Confinement_Pressure)
*format "%10.5e"
InitialTime = *GenData(Initial_Time)
*format "%10.5e"
IncreasingTemporaily = *GenData(Increasing_Temporaily)

# Rotation Options

RotationOption= "*GenData(Rotation_Option)"
TrihedronOption= "*GenData(Trihedron_Option)"
RotationalSpringOption= "*GenData(Rotation_Spring)"
RotaDampId = "*GenData(Rota_Damp_Id)"

#POSTPROCES

print_velocity                   = "*GenData(VELOCITY)"
print_displacement               = "*GenData(DISPLACEMENT)"
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


mass_elements=*tcl(DEM::Get_Mass_Elements)

# Declare Python Variables
