
# DEM General Options
Dimension                        = 3
PeriodicDomainOption             = "OFF"
BoundingBoxOption                = "ON"
AutomaticBoundingBoxOption       = "ON"
BoundingBoxEnlargementFactor     = 1.1
BoundingBoxStartTime             = 0.0
BoundingBoxStopTime              = 1000.0
BoundingBoxMaxX                  =  1.00000e+01
BoundingBoxMaxY                  =  1.00000e+01
BoundingBoxMaxZ                  =  1.00000e+01
BoundingBoxMinX                  = -1.00000e+01
BoundingBoxMinY                  = -1.00000e+01
BoundingBoxMinZ                  = -1.00000e+01

dem_inlet_option                 = 0
GravityX                         = 0.0
GravityY                         = 0.0
GravityZ                         = -0.0

VelocityTrapOption               = 0
RotationOption                   = "ON"
CleanIndentationsOption          = "OFF"
DeltaOption                      = "Absolute"
SearchTolerance                  = 0.01
CoordinationNumber               = 10
AmplifiedSearchRadiusExtension   = 1e-4
ModelDataInfo                    = "OFF"
VirtualMassCoefficient           = 1.0
RollingFrictionOption            = "OFF"
DontSearchUntilFailure           = "No"
ContactMeshOption                = "ON"
OutputFileType                   = "Binary"
Multifile                        = "multiple_files"
MaxAmplificationRatioOfSearchRadius = 1000

# Solution Strategy

IntegrationScheme                = "Symplectic_Euler"
AutomaticTimestep                = "OFF"
DeltaTimeSafetyFactor            = 1.0
MaxTimeStep                      = 5e-7
FinalTime                        = 0.05
ControlTime                      = 4.0
NeighbourSearchFrequency         = 50

# Material Test
TestType                         = "None"
ConfinementPressure              = 0.0
LoadingVelocityTop               = 0.0
LoadingVelocityBot               = 0.0
StressStrainOption               = 1
MeshType                         = "Current"
MeshPath                         = "0"
SpecimenLength                   = 0.30
SpecimenDiameter                 = 0.15
MeasuringSurface                 = 0.01767145867644375

ElementType                      = "SphericContPartDEMElement3D"


# PostProcess Results

GraphExportFreq                  = 1e-3
VelTrapGraphExportFreq           = 1e0
OutputTimeStep                   = 0.5e-3
Granulometry                     = "No"
PostDisplacement                 = 1
PostVelocity                     = 1
PostElasticForces                = 1
PostContactForces                = 1
PostRigidElementForces           = 0
PostTangentialElasticForces      = 0
PostTotalForces                  = 1
PostShearStress                  = 0
PostNonDimensionalVolumeWear     = 0
PostNodalArea                    = 0
PostStressStrainOption           = 0
PostTemperature                  = 0
PostHeatFlux                     = 0
PostSkinSphere                   = 0
PostPoissonRatio                 = 0
PostRHS                          = 0
PostDampForces                   = 0
PostAppliedForces                = 0
PostRadius                       = 1
PostGroupId                      = 0
PostExportId                     = 0
PostAngularVelocity              = 1
PostParticleMoment               = 0
PostEulerAngles                  = 0
PostContactSigma                 = 0
PostContactTau                   = 0
PostLocalContactForce            = 0
PostFailureCriterionState        = 1
PostContactFailureId             = 1
PostMeanContactArea              = 0
PostBoundingBox                  = 0

PostPressure                     = 0
# FROM CND:
PredefinedSkinOption             = "OFF"
MeanRadius                       = 0.0001

# Declare Python Variables
problem_name="ucs_benchmark"
kratos_path="D:\Kratos"
