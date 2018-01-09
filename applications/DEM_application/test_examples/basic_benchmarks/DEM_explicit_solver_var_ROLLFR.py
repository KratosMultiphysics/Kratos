
# DEM General Options

Dimension                        = 3
BoundingBoxOption                = "ON"
BoundingBoxEnlargementFactor     = 1.1
AutomaticBoundingBoxOption       = "OFF"
BoundingBoxEnlargementFactor     = 1.0
BoundingBoxStartTime             = 0.0
BoundingBoxStopTime              = 0.7001
BoundingBoxMaxX                  = 0.008
BoundingBoxMaxY                  = 0.008
BoundingBoxMaxZ                  = 0.003
BoundingBoxMinX                  = -0.003
BoundingBoxMinY                  = 0
BoundingBoxMinZ                  = -0.003

dem_inlet_option                 = 0
GravityX                         = 0.0
GravityY                         = -9.81
GravityZ                         = 0.0

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
MaxTimeStep                      = 5e-7
FinalTime                        = 0.2001
ControlTime                      = 100.0
NeighbourSearchFrequency         = 1
PeriodicDomainOption             = 0
# Material Test

TestType                         = "None"
ConfinementPressure              = 0.0
LoadingVelocityTop               = 0.0
LoadingVelocityBot               = 0.0
MeshType                         = "Current"
MeshPath                         = "0"
SpecimenLength                   = 0.30
SpecimenDiameter                 = 0.15
MeasuringSurface                 = 0.01767145867644375

ElementType                      = "SphericPartDEMElement3D"
# PostProcess

# PostProcess Results

GraphExportFreq                  = 1e-3
VelTrapGraphExportFreq           = 1e-3
OutputTimeStep                   = 0.01
PostDisplacement                 = "1"
PostVelocity                     = "1"
PostElasticForces                = "0"
PostContactForces                = "0"
PostRigidElementForces           = "0"
PostTangentialElasticForces      = "0"
PostPressure                     = "0"
PostTotalForces                  = "0"
StressStrainOption               = "0"
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
PostParticleMoment               = 1
PostEulerAngles                  = 1
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

problem_name="rolling_benchmark"
kratos_path="D:\Kratos"
