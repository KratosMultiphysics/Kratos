
# DEM General Options

Dimension                        = 3
BoundingBoxOption                = "OFF"
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
RollingFrictionOption            = "OFF"
DontSearchUntilFailure           = "OFF"
ContactMeshOption                = "OFF"
OutputFileType                   = "Binary"
Multifile                        = "multiple_files"
HorizontalFixVel                 = "ON"

# Solution Strategy

IntegrationScheme                = "Taylor_Scheme"
AutomaticTimestep                = "OFF"
DeltaTimeSafetyFactor            = 1.0
MaxTimeStep                      = 1e-5
FinalTime                        = 2.0
ControlTime                      = 100.0
NeighbourSearchFrequency         = 1
PeriodicDomainOption             = 0
# Material Test

TestType                         = "None"
ConfinementPressure              = 0.0
LoadingVelocityTop               = 0.0
LoadingVelocityBot               = 0.0
StressStrainOption               = "OFF"
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
PostDisplacement                 = "0"
PostVelocity                     = "0"
PostElasticForces                = "0"
PostContactForces                = "0"
PostRigidElementForces           = "0"
PostTangentialElasticForces      = "0"
PostPressure                     = "0"
PostTotalForces                  = "0"
PostStressStrainOption           = "0"
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
PostAngularVelocity              = 0
PostParticleMoment               = 0
PostEulerAngles                  = 0
PostContactSigma                 = 0
PostContactTau                   = 0
PostLocalContactForce            = 0
PostFailureCriterionState        = 0
PostContactFailureId             = 0
PostMeanContactArea              = 0
PostBoundingBox                  = 0

# FROM CND:

PredefinedSkinOption             = "OFF"
MeanRadius                       = 0.0001

# Declare Python Variables

problem_name="benchmark13"
kratos_path="D:\Kratos"
