# Model Type

ContinuumOption                  = "OFF"
RotationOption                   = "OFF"
HomogeneousMaterialOption        = "OFF"
ElementType                      = "SphericParticle3D"

# Meshing Settings

CleanIndentationsOption          = "OFF"

# General Settings

final_time                       = 1.00000e+00
gravity_x                        = 0.00000e+00
gravity_y                        = -9.81000e+00
gravity_z                        = 0.00000e+00
BoundingBoxOption                = "OFF"
bounding_box_enlargement_factor  = 2.00000e+00
domain_size                      = 3
OutputFileType                   = "Binary"
Multifile                        = "single_file"
PrintNeighbourLists              = "OFF"
ModelDataInfo                    = "OFF"

# Special features

VirtualMassOption                = "OFF"
VirtualMassCoefficient           = 0.00000e+00
MagicFactor                      = 1.00000e+00
DeltaOption                      = "OFF"
search_radius_extension          = 1.00000e-02
ampl_search_radius_extension     = 1.10000e+00
FixVelocities                    = "OFF"
TimePercentageFixVelocities      = 1.00000e+01
TrihedronOption                  = "OFF"
LimitSurfaceOption               = "OFF"
surface_normal_dir_x             = 0.00000e+00
surface_normal_dir_y             = 1.00000e+00
surface_normal_dir_z             = 0.00000e+00
surface_point_coor_x             = 0.00000e+00
surface_point_coor_y             = 0.00000e+00
surface_point_coor_z             = 0.00000e+00
surface_friction_angle           = 45

# Time Discretization Settings

Integration_Scheme               = "forward_euler"
search_step                      = 1.00000e+00
AutoReductionOfTimeStepOption    = "ON"
dt_safety_factor                 = 1.00000e+00
max_time_step                    = 1.00000e-04
output_dt                        = 1.00000e-02
control_time                     = 2.00000e+01

# Material Model

NormalForceCalculation           = "Linear"
NormalDampId                     = "ViscDamp"
TangentialDampId                 = "NoDamp"
FailureCriterionOption           = "Uncoupled"
TauZero                          = 5.26000e+00
SigmaMax                         = 3.32000e+01
SigmaMin                         = 3.32000e+00
InternalFricc                    = 3.50000e+01
NonLinearOption                  = "OFF"
C1                               = 4.00000e-01
N1                               = 2.00000e+01
C2                               = 1.50000e+00
N2                               = 1.50000e+01
RotationalSpringOption           = "OFF"
RotaDampId                       = "NoDamp"

# Global Material Parameters

GeneralDensity                   = "1.00000e+02"
GeneralYoungModulus              = "2.10000e+07"
GeneralPoissonRatio              = "5.00000e-01"
GeneralCohesion                  = "4.16000e+06"
GeneralRollingFriction           = "0.00000e+00"
GeneralTension                   = "2.01000e+06"
GeneralRotaDampRatio             = "5.00000e-01"
GeneralStaticFrictionCoef        = "0.00000e+00"
GeneralDynamicFrictionCoef       = "0.00000e+00"
GeneralRestitutionCoef           = "5.00000e-01"
GeneralColour                    = "1.00000e+00"
GlobalVariablesOption            = "OFF"
global_kn                        = 3.00000e+03
global_kt                        = 1.00000e+03
global_kr                        = 1.00000e+03
global_rn                        = 1.00000e+03
global_rt                        = 1.00000e+03
global_rr                        = 1.00000e+03
global_fri_ang                   = 4.00000e+01

# Continuum Options

StressStrainOperations           = "OFF"
ContactMeshOption                = "OFF"
ConcreteTestOption               = "OFF"
RealTimeGraph                    = "OFF"
TriaxialOption                   = "OFF"
ConfinementPressure              = 0.00000e+00
InitialTime                      = 0.00000e+00
IncreasingTemporaily             = 1.50000e+01

#POSTPROCES

print_velocity                   = "1"
print_displacement               = "1"
print_radial_displacement        = "0"
print_rhs                        = "0"
print_total_forces               = "0"
print_damp_forces                = "0"
print_applied_forces             = "0"
print_radius                     = "1"
print_particle_cohesion          = "0"
print_particle_tension           = "0"
print_group_id                   = "0"
print_export_id                  = "1"
print_export_particle_failure_id = "0"
print_export_skin_sphere         = "0"
print_local_contact_force_low    = "0"
print_local_contact_force_high   = "0"
print_failure_criterion_state    = "0"
print_contact_failure            = "0"
print_contact_tau                = "0"
print_contact_sigma              = "0"
print_angular_velocity           = "0"
print_particle_moment            = "0"
print_euler_angles               = "0"
print_representative_volume      = "0"
print_mean_contact_area          = "0"
print_stress_tensor              = "0"

#FROM CND:

predefined_skin_option           = "OFF"

mass_elements                    = 2.26195e-01

# Declare Python Variables
problem_name = 'two_balls_normal_damp'
problem_path = '/home/cimne/kratos/applications/DEM_application/test_examples/two_balls_normal_damp.gid'
kratos_path = '/home/cimne/kratos'

