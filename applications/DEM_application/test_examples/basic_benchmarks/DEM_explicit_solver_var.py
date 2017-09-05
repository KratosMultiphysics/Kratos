
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
GravityZ                         = 0.0

VelocityTrapOption               = 0
RotationOption                   = "ON"

Dempack                          = "OFF"
CleanIndentationsOption          = "ON"

RemoveBallsInEmbeddedOption      = 1

DeltaOption                      = "Absolute"
SearchTolerance                  = 0.0
CoordinationNumber               = 10
AmplifiedSearchRadiusExtension   = 1
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
MaxTimeStep                      = 1
FinalTime                        = 1
ControlTime                      = 100
NeighbourSearchFrequency         = 1
PeriodicDomainOption             = 0

# Material Test

TestType                         = "None"
ConfinementPressure              = 0.0
LoadingVelocityTop               = 0.0
LoadingVelocityBot               = 0.0
FemPlates                        = "OFF"
StressStrainOption               = "OFF"
MeshType                         = "Current"
MeshPath                         = "0"
SpecimenLength                   = 1
SpecimenDiameter                 = 1
MeasuringSurface                 = 1

ElementType                      = "SphericPartDEMElement3D"

# PostProcess Results

GraphExportFreq                  = 100
VelTrapGraphExportFreq           = 100
OutputTimeStep                   = 100
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

problem_name="DEFAULT"
kratos_path="D:\Kratos"
