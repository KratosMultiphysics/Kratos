
# DEM General Options

Dimension                        = 3
drag_modifier_type               = 3
project_from_particles_option    = 0
consider_lift_force_option       = 0
BoundingBoxOption                = "ON"
BoundingBoxEnlargementFactor     = 1.1
AutomaticBoundingBoxOption       = "OFF"
BoundingBoxEnlargementFactor     = 1.0
BoundingBoxMaxX                  = 3.35
BoundingBoxMaxY                  = 1.18
BoundingBoxMaxZ                  = 6
BoundingBoxMinX                  = -4.1
BoundingBoxMinY                  = -1.88
BoundingBoxMinZ                  = -0.1

dem_inlet_option                 = 1
GravityX                         = 0.0
GravityY                         = 0.0
GravityZ                         = -9.81

VelocityTrapOption               = 0
RotationOption                   = "ON"

Dempack                          = "OFF"
CleanIndentationsOption          = "ON"

RemoveBallsInEmbeddedOption      = 1

DeltaOption                      = "Absolute"
SearchTolerance                  = 0.0
CoordinationNumber               = 10
AmplifiedSearchRadiusExtension   = 1.10000e+00
ModelDataInfo                    = "OFF"
VirtualMassCoefficient           = 1.0
RollingFrictionOption            = "ON"
DontSearchUntilFailure           = "OFF"
ContactMeshOption                = "OFF"
OutputFileType                   = "Binary"
Multifile                        = "multiple_files"
HorizontalFixVel                 = "ON"

# Solution Strategy

IntegrationScheme                = "Forward_Euler"
AutomaticTimestep                = "OFF"
DeltaTimeSafetyFactor            = 1.0
MaxTimeStep                      = 5e-5
FinalTime                        = 15.0
ControlTime                      = 2.0
NeighbourSearchFrequency         = 1

# Constitutive Parameters

MaterialModel                    = "Linear"
G1                               = 0.0
G2                               = 0.0
G3                               = 0.0
MaxDef                           = 0.0
FailureCriterionType             = "Uncoupled"
AreaFactor                       = 1.00000e+00
LocalContactDamping              = "Normal"
LocalDampingFactor               = 0.9
GlobalForceReduction             = 0.2

# Material Test

TestType                         = "None"
ConfinementPressure              = 0.0
LoadingVelocityTop               = 0.0
LoadingVelocityBot               = 0.0
FemPlates                        = "OFF"
StressStrainOption               = "OFF"
MeshType                         = "Current"
MeshPath                         = "0"
SpecimenLength                   = 0.30
SpecimenDiameter                 = 0.15
MeasuringSurface                 = 0.01767145867644375

ElementType                      = "SphericPartDEMElement3D"

# PostProcess Results

GraphExportFreq                  = 1e-6
VelTrapGraphExportFreq           = 1e-3
OutputTimeStep                   = 0.05
PostDisplacement                 = "1"
PostVelocity                     = "1"
PostElasticForces                = "0"
PostContactForces                = "0"
PostRigidElementForces           = "0"
PostTangentialElasticForces      = "0"
PostPressure                     = "0"
PostTotalForces                  = "0"
PostShearStress                  = "0"
PostNonDimensionalVolumeWear     = "0"
PostNodalArea                    = "0"
PostRHS                          = "0"
PostDampForces                   = "0"
PostAppliedForces                = "0"
PostRadius                       = "0"
PostGroupId                      = "0"
PostExportId                     = "0"
PostExportSkinSphere             = "0"
PostAngularVelocity              = 1
PostParticleMoment               = 0
PostEulerAngles                  = 1
PostContactSigma                 = 0
PostContactTau                   = 0
PostLocalContactForce            = 0
PostFailureCriterionState        = 0
PostContactFailureId             = 0
PostMeanContactArea              = 0

# FROM CND:

PredefinedSkinOption             = "OFF"
MeanRadius                       = 0.0001

# Declare Python Variables

problem_name="D-Dempack_validation"
kratos_path="D:\Kratos"
