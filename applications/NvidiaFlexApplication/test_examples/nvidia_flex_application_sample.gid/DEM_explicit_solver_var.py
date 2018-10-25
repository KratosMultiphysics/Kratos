
# DEM General Options
Dimension                        = 3
PeriodicDomainOption             = "OFF"
BoundingBoxOption                = "ON"
AutomaticBoundingBoxOption       = "OFF"
BoundingBoxEnlargementFactor     = 1.0
BoundingBoxStartTime             = 0.0
BoundingBoxStopTime              = 1000.0
BoundingBoxMaxX                  = 0.0275
BoundingBoxMaxY                  = 0.0275
BoundingBoxMaxZ                  = 0.057
BoundingBoxMinX                  = -0.0275
BoundingBoxMinY                  = -0.0275
BoundingBoxMinZ                  = 0.042

dem_inlet_option                 = 0
GravityX                         = 0.0
GravityY                         = 0.0
GravityZ                         = -9.81

EnergyCalculationOption          = 0
VelocityTrapOption               = 0
RotationOption                   = "OFF"
CleanIndentationsOption          = "OFF"
RemoveBallsInEmbeddedOption      = 1

DeltaOption                      = "Absolute"
SearchTolerance                  = 0.0
AmplifiedSearchRadiusExtension   = 0.0
ModelDataInfo                    = "OFF"
VirtualMassCoefficient           = 1.0
RollingFrictionOption            = "OFF"
ContactMeshOption                = "OFF"
OutputFileType                   = "Binary"
Multifile                        = "multiple_files"
ElementType                      = "SphericPartDEMElement3D"

# Solution Strategy
TranslationalIntegrationScheme       = "Symplectic_Euler"
RotationalIntegrationScheme          = "Direct_Integration"
AutomaticTimestep                = "OFF"
DeltaTimeSafetyFactor            = 1.0
MaxTimeStep                      = 1.5e-5
FinalTime                        = 2.5
ControlTime                      = 4.0
NeighbourSearchFrequency         = 50

# PostProcess Results
GraphExportFreq                  = 1e-3
VelTrapGraphExportFreq           = 1e-3
OutputTimeStep                   = 0.02
PostBoundingBox                  = 0
PostDisplacement                 = 0
PostVelocity                     = 0
# DEM only Results
PostTotalForces                  = 0
PostRigidElementForces           = 0
PostRadius                       = 0
PostAngularVelocity              = 0
PostParticleMoment               = 0
PostEulerAngles                  = 0
PostRollingResistanceMoment      = 0
# FEM only Results
PostElasticForces                = 0
PostContactForces                = 0
PostTangentialElasticForces      = 0
PostShearStress                  = 0
PostPressure                     = 0
# FEM_wear only Results
PostNonDimensionalVolumeWear     = 0
PostNodalArea                    = 0
# Results on bond elements
# Under revision
PostRHS                          = 0
PostDampForces                   = 0
PostAppliedForces                = 0
PostGroupId                      = 0
PostExportId                     = 0

#
problem_name="nvidia_flex_application_sample"
