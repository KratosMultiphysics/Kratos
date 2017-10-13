
# DEM General Options
Dimension                        = 3
PeriodicDomainOption             = "OFF"
BoundingBoxOption                = "ON"
BoundingBoxEnlargementFactor     = 1.0
AutomaticBoundingBoxOption       = "OFF"
BoundingBoxMaxX                  =  1.00000e+01
BoundingBoxMaxY                  =  1.00000e+01
BoundingBoxMaxZ                  =  1.00000e+01
BoundingBoxMinX                  = -1.00000e+01
BoundingBoxMinY                  = -1.00000e+01
BoundingBoxMinZ                  = -1.00000e+01

dem_inlet_option                 = 0
GravityX                         = 0.0
GravityY                         = 0.0
GravityZ                         = 0.0

VelocityTrapOption               = 0
RotationOption                   = "ON"
CleanIndentationsOption          = "OFF"
DeltaOption                      = "Absolute"
SearchTolerance                  = 0.001
CoordinationNumber               = 10
AmplifiedSearchRadiusExtension   = 1.1
ModelDataInfo                    = "OFF"
VirtualMassCoefficient           = 1.0
RollingFrictionOption            = "OFF"
PoissonEffectOption              = "OFF"
ShearStrainParallelToBondOption  = "OFF"
DontSearchUntilFailure           = "OFF"
ContactMeshOption                = "ON"
OutputFileType                   = "Binary"
Multifile                        = "multiple_files"

# Solution Strategy

IntegrationScheme                = "Forward_Euler" #"Verlet_Velocity"   #"Symplectic_Euler" Forward_Euler
AutomaticTimestep                = "OFF"
DeltaTimeSafetyFactor            = 1.0
MaxTimeStep                      = 1e-5
FinalTime                        = 0.01
ControlTime                      = 100.0
NeighbourSearchFrequency         = 1000

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
MaxAmplificationRatioOfSearchRadius = 1000

ElementType                      = "SphericContPartDEMElement3D"

# PostProcess Results

GraphExportFreq                  = 1e-5
VelTrapGraphExportFreq           = 1e0
OutputTimeStep                   = 1e-3
Granulometry                     = "No"
PostDisplacement                 = 1
PostVelocity                     = 1
PostElasticForces                = 1
PostContactForces                = 0
PostRigidElementForces           = 0
PostTangentialElasticForces      = 0
PostTotalForces                  = 1
PostShearStress                  = 1
PostNonDimensionalVolumeWear     = 0
PostNodalArea                    = 0
PostStressStrainOption           = 1
PostTemperature                  = 0
PostHeatFlux                     = 0
PostSkinSphere                   = 1
PostPoissonRatio                 = 0
PostRHS                          = 0
PostDampForces                   = 0
PostAppliedForces                = 0
PostRadius                       = 0
PostGroupId                      = 0
PostExportId                     = 0
PostAngularVelocity              = 1
PostParticleMoment               = 1
PostEulerAngles                  = 0
PostContactSigma                 = 1
PostContactTau                   = 1
PostLocalContactForce            = 1
PostFailureCriterionState        = 0
PostContactFailureId             = 1
PostMeanContactArea              = 0
PostBoundingBox                  = 0
PostPressure                     = 0

#
problem_name="benchmark20"
