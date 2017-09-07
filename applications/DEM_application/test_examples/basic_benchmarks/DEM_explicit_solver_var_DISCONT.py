
# DEM General Options

Dimension                        = 3
drag_modifier_type               = 3
project_from_particles_option    = 0
consider_lift_force_option       = 0
BoundingBoxOption                = "ON"
BoundingBoxEnlargementFactor     = 1.1
AutomaticBoundingBoxOption       = "OFF"
BoundingBoxEnlargementFactor     = 1.0
BoundingBoxMaxX                  = 1e3
BoundingBoxMaxY                  = 1e3
BoundingBoxMaxZ                  = 1e3
BoundingBoxMinX                  = -1e3
BoundingBoxMinY                  = -1e3
BoundingBoxMinZ                  = -1e3

dem_inlet_option                 = 0
GravityX                         = 0.0
GravityY                         = 0.0
GravityZ                         = 0.0 #-9.81

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
RollingFrictionOption            = "OFF"
DontSearchUntilFailure           = "OFF"
ContactMeshOption                = "OFF"
OutputFileType                   = "Binary"
Multifile                        = "multiple_files"
HorizontalFixVel                 = "ON"

# Solution Strategy

IntegrationScheme                = "Forward_Euler"
AutomaticTimestep                = "OFF"
DeltaTimeSafetyFactor            = 1.0
MaxTimeStep                      = 6.4e-8
FinalTime                        = 0.0005
ControlTime                      = 100
NeighbourSearchFrequency         = 1
PeriodicDomainOption             = 0

# Constitutive Parameters

MaterialModel                    = "Hertz"
G1                               = 0.0
G2                               = 0.0
G3                               = 0.0
MaxDef                           = 0.0
FailureCriterionType             = "Uncoupled"
AreaFactor                       = 1.00000e+00
LocalContactDamping              = "Normal"
LocalDampingFactor               = 1.0 #0.9
GlobalForceReduction             = 0.0 #0.2

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

GraphExportFreq                  = 1e-5
VelTrapGraphExportFreq           = 1e-3
OutputTimeStep                   = 0.05
PostDisplacement                 = "1"
PostVelocity                     = "1"
PostElasticForces                = "1"
PostContactForces                = "1"
PostRigidElementForces           = "1"
PostTangentialElasticForces      = "0"
PostPressure                     = "1"
PostTotalForces                  = "1"
PostShearStress                  = "0"
PostNonDimensionalVolumeWear     = "0"
PostNodalArea                    = "0"
PostRHS                          = "1"
PostDampForces                   = "0"
PostAppliedForces                = "1"
PostRadius                       = "1"
PostGroupId                      = "0"
PostExportId                     = "0"
PostExportSkinSphere             = "0"
PostAngularVelocity              = 1
PostParticleMoment               = 0
PostEulerAngles                  = 0
PostContactSigma                 = 0
PostContactTau                   = 0
PostLocalContactForce            = 0
PostFailureCriterionState        = 0
PostContactFailureId             = 0
PostMeanContactArea              = 0
PostStressStrainOption           = 0

# FROM CND:

PredefinedSkinOption             = "OFF"
MeanRadius                       = 0.0001

# Declare Python Variables

problem_name="benchmark5"
kratos_path="D:\Kratos"
