# General Settings

BoundingBoxOption= "OFF"
bounding_box_enlargement_factor = 2.00000e+00
domain_size = 3
OutputFileType="Binary"
Multifile="single_file"

# Calculation Settings

Integration_Scheme = "forward_euler"
NormalForceCalculation = "Linear"
DampId = "NoDamp"
gravity_x = 0.00000e+00
gravity_y = -9.81000e+00
gravity_z = 0.00000e+00

#Time Settings

search_step = 1.00000e+00
final_time = 3.00000e+00
CriticalTimeOption = "ON"
dt_safety_factor = 1.00000e+00
max_time_step = 1.00000e-04
output_dt = 1.00000e-02
control_time = 2.00000e+01

#Special features

VirtualMassOption = "OFF"
VirtualMassCoefficient = 0.00000e+00
MagicFactor = 1.00000e+00
DeltaOption = "OFF"
search_radius_extension = 1.00000e-02
PrintNeighbourLists= "OFF"
GlobalVariablesOption= "OFF"
global_kn = 3.00000e+03
global_kt = 1.00000e+03
global_kr = 1.00000e+03
global_rn = 1.00000e+03
global_rt = 1.00000e+03
global_rr = 1.00000e+03
global_fri_ang = 4.00000e+01
ModelDataInfo = "OFF"

#Continuum Options

ContinuumOption= "OFF"
ContactMeshOption = "OFF"
FailureCriterionOption = "Uncoupled"
TauZero = 5.26000e+00
SigmaMax = 3.32000e+01
SigmaMin = 3.32000e+00
InternalFricc = 3.50000e+01

#Concrete Test Options

ConcreteTestOption = "OFF"
RealTimeGraph = "OFF"
TriaxialOption = "OFF"
ConfinementPressure = 0.00000e+00
InitialTime = 0.00000e+00
IncreasingTemporaily = 1.50000e+01

# Rotation Options

RotationOption= "OFF"
TrihedronOption= "OFF"
RotationalSpringOption= "OFF"
RotaDampId = "NoDamp"

#POSTPROCES

print_velocity                   = "1"
print_displacement               = "1"
print_rhs                        = "1"
print_total_forces               = "0"
print_damp_forces                = "0"
print_applied_forces             = "0"
print_radius                     = "0"
print_particle_cohesion          = "0"
print_particle_tension           = "0"
print_group_id                   = "0"
print_export_id                  = "1"
print_export_particle_failure_id = "0"
print_export_skin_sphere         = "1"
print_local_contact_force_low    = "0"
print_local_contact_force_high   = "0"
print_failure_criterion_state    = "1"
print_contact_failure            = "1"
print_contact_tau                = "1"
print_contact_sigma              = "1"
print_angular_velocity           = "0"
print_particle_moment            = "0"
print_euler_angles               = "0"


mass_elements=1.30900e-01

# Declare Python Variables
problem_name = 'two_balls_no_damp'
problem_path = '/home/CIMNE/kratos/applications/DEM_application/test_examples/two_balls_no_damp.gid'
kratos_path = '/home/CIMNE/kratos'

