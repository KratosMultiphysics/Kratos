
# DEM General Options
Dimension                        = 3
PeriodicDomainOption             = "OFF"
BoundingBoxOption                = "ON"
AutomaticBoundingBoxOption       = "OFF"
BoundingBoxEnlargementFactor     = 1.0
BoundingBoxStartTime             = 0.0
BoundingBoxStopTime              = 1000.0
BoundingBoxMaxX                  = 1000
BoundingBoxMaxY                  = 1000
BoundingBoxMaxZ                  = 1000
BoundingBoxMinX                  = -1000
BoundingBoxMinY                  = -1000
BoundingBoxMinZ                  = -1000

dem_inlet_option                 = 0
GravityX                         = 0.0
GravityY                         = 0.0
GravityZ                         = -9.81

VelocityTrapOption               = 0
RotationOption                   = "ON"
CleanIndentationsOption          = "OFF"
RemoveBallsInEmbeddedOption      = 1

DeltaOption                      = "Absolute"
SearchTolerance                  = 0.0001
CoordinationNumber               = 10
AmplifiedSearchRadiusExtension   = 0.0
ModelDataInfo                    = "OFF"
VirtualMassCoefficient           = 1.0
RollingFrictionOption            = "OFF"
DontSearchUntilFailure           = "OFF"
ContactMeshOption                = "OFF"
OutputFileType                   = "Binary"
Multifile                        = "multiple_files"

# Solution Strategy

IntegrationScheme                = "Hybrid_Bashforth"
AutomaticTimestep                = "OFF"
DeltaTimeSafetyFactor            = 1.0
MaxTimeStep                      = 0.01
FinalTime                        = 1.0
ControlTime                      = 4.0
NeighbourSearchFrequency         = 1
TestType = "None"

ElementType                      = "SwimmingDEMElement"

# Swimming DEM-specific section begins
#-------------------------------------
coupling_level_type                    = 1
time_averaging_type                    = 0
pick_individual_forces_option          = 0
include_faxen_terms_option             = 0 # (relevant if the Maxey Riley equation is used)
gradient_calculation_type              = 1 # (Not calculated (0), volume-weighed average(1), Superconvergent recovery(2))
laplacian_calculation_type             = 2 # (Not calculated (0), Finite element projection (1), Superconvergent recovery(2))
buoyancy_force_type                    = 2 # null buoyancy (0), compute buoyancy (1)  if drag_force_type is 2 buoyancy is always parallel to gravity
drag_force_type                        = 10 # null drag (0), Stokes (1), Weatherford (2), Ganser (3), Ishii (4), Newtonian Regime (5)
virtual_mass_force_type                = 10 # null virtual mass force (0)
lift_force_type                        = 0 # null lift force (0), Saffman (1)
magnus_force_type                      = 0 # null magnus force (0), Rubinow and Keller (1), Oesterle and Bui Dihn (2)
hydro_torque_type                      = 0 # null hydrodynamic torque (0), Dennis (1)
drag_modifier_type                     = 0
viscosity_modification_type            = 0

# Parameters not yet settable from interface

coupling_weighing_type                 = 2 # {fluid_to_DEM, DEM_to_fluid, fluid_fraction} = {lin, lin, imposed} (-1), {lin, const, const} (0), {lin, lin, const} (1), {lin, lin, lin} (2), averaging method (3)
fluid_model_type                       = 1 # untouched, velocity incremented by 1/fluid_fraction (0), modified mass conservation only (1)
coupling_scheme_type                   = "UpdatedFluid" # "UpdatedFluid", "UpdatedDEM"
print_particles_results_option         = 0
add_each_hydro_force_option            = 1 # add each of the hydrodynamic forces (drag, lift and virtual mass)
project_at_every_substep_option        = 1
velocity_trap_option                   = 0
inlet_option                           = 1
manually_imposed_drag_law_option       = 0
stationary_problem_option              = 0 # stationary, stop calculating the fluid after it reaches the stationary state
flow_in_porous_medium_option           = 0 # the porosity is an imposed field
flow_in_porous_DEM_medium_option       = 0 # the DEM part is kept static
embedded_option                        = 0 # the embedded domain tools are to be used
make_results_directories_option        = 1 # results are written into a folder (../results) inside the problem folder
body_force_on_fluid_option             = 1
print_debug_info_option                = 0 # print a summary of global physical measures
print_particles_results_cycle          = 1 # number of 'ticks' per printing cycle
debug_tool_cycle                       = 10 # number of 'ticks' per debug computations cycle
similarity_transformation_type         = 0 # no transformation (0), Tsuji (1)
dem_inlet_element_type                 = "SphericSwimmingParticle3D"  # "SphericParticle3D", "SphericSwimmingParticle3D"
drag_modifier_type                     = 2 # Hayder (2), Chien (3) # problemtype option
drag_porosity_correction_type          = 0 # No correction (0), Richardson and Zaki (1)
interaction_start_time                 = 0
min_fluid_fraction                     = 0.2
initial_drag_force                     = 0.0   # problemtype option
drag_law_slope                         = 0.0   # problemtype option
power_law_tol                          = 0.0
model_over_real_diameter_factor        = 1.0 # not active if similarity_transformation_type = 0
max_pressure_variation_rate_tol        = 1e-3 # for stationary problems, criterion to stop the fluid calculations
time_steps_per_stationarity_step       = 15 # number of fluid time steps between consecutive assessment of stationarity steps
meso_scale_length                      = 0.2 # the radius of the support of the averaging function for homogenization (<=0 for automatic calculation)
shape_factor                           = 0.5 # the density function's maximum over its support's radius (only relevant if coupling_weighing_type == 3)
#-------------------------------------
# Swimming DEM-specific section ends



# PostProcess Results

GraphExportFreq                  = 1e-3
VelTrapGraphExportFreq           = 1e-3
OutputTimeStep                   = 0.5
PostDisplacement                 = 1
PostVelocity                     = 1
PostElasticForces                = 0
PostContactForces                = 0
PostRigidElementForces           = 0
PostTangentialElasticForces      = 0
PostTotalForces                  = 0
PostShearStress                  = 0
PostNonDimensionalVolumeWear     = 0
PostNodalArea                    = 0
PostRHS                          = 0
PostDampForces                   = 0
PostAppliedForces                = 0
PostRadius                       = 1
PostGroupId                      = 0
PostExportId                     = 0
PostAngularVelocity              = 0
PostParticleMoment               = 0
PostEulerAngles                  = 0
PostBoundingBox                  = 0

PostPressure                     = 0
# Swimming DEM-specific section begins
#-------------------------------------
PostFluidPressure                          = 0
print_REYNOLDS_NUMBER_option               = 0
print_PRESSURE_GRAD_PROJECTED_option       = 0
print_FLUID_VEL_PROJECTED_option           = 1
print_FLUID_ACCEL_PROJECTED_option         = 1
print_BUOYANCY_option                      = 1
print_DRAG_FORCE_option                    = 1
print_VIRTUAL_MASS_FORCE_option            = 1
print_BASSET_FORCE_option                  = 1
print_LIFT_FORCE_option                    = 1
print_FLUID_VEL_PROJECTED_RATE_option      = 1
print_FLUID_VISCOSITY_PROJECTED_option     = 0
print_FLUID_FRACTION_PROJECTED_option      = 0
print_FLUID_VEL_LAPL_PROJECTED_option      = 0
print_FLUID_VEL_LAPL_RATE_PROJECTED_option = 0
print_HYDRODYNAMIC_FORCE_option            = 1
print_HYDRODYNAMIC_MOMENT_option           = 0
print_MESH_VELOCITY1_option                = 0
print_BODY_FORCE_option                    = 0
print_FLUID_FRACTION_option                = 0
print_FLUID_FRACTION_GRADIENT_option       = 0
print_HYDRODYNAMIC_REACTION_option         = 0
print_PRESSURE_option                      = 0
print_PRESSURE_GRADIENT_option             = 0
print_DISPERSE_FRACTION_option             = 0
print_MEAN_HYDRODYNAMIC_REACTION_option    = 0
print_VELOCITY_LAPLACIAN_option            = 0
print_VELOCITY_LAPLACIAN_RATE_option       = 0

#-------------------------------------
# Swimming DEM-specific section ends

# FROM CND:
PredefinedSkinOption             = "OFF"
MeanRadius                       = 0.0001

# Declare Python Variables
problem_name="Candelier"
kratos_path="D:\kratos"
