
# DEM General Options
Dimension                        = 3
PeriodicDomainOption             = "OFF"
BoundingBoxOption                = "OFF"
AutomaticBoundingBoxOption       = "OFF"
BoundingBoxEnlargementFactor     = 1.0
BoundingBoxStartTime             = 0.0
BoundingBoxStopTime              = 1000.0
BoundingBoxMaxX                  = 10
BoundingBoxMaxY                  = 10
BoundingBoxMaxZ                  = 10
BoundingBoxMinX                  = -10
BoundingBoxMinY                  = -10
BoundingBoxMinZ                  = -10

dem_inlet_option                 = 0
GravityX                         = 0.0
GravityY                         = 0.0
GravityZ                         = -9.81

EnergyCalculationOption          = 1
PotentialEnergyReferencePointX   = 0.0
PotentialEnergyReferencePointY   = 0.0
PotentialEnergyReferencePointZ   = 0.0
VelocityTrapOption               = 0
RotationOption                   = "ON"
CleanIndentationsOption          = "OFF"
RemoveBallsInEmbeddedOption      = 1

DeltaOption                      = "Absolute"
SearchTolerance                  = 0.0
AmplifiedSearchRadiusExtension   = 0.0
ModelDataInfo                    = "OFF"
VirtualMassCoefficient           = 1.0
RollingFrictionOption            = "OFF"
ContactMeshOption                = "OFF"
OutputFileType                   = "Ascii"
Multifile                        = "multiple_files"
ElementType                      = "SphericPartDEMElement3D"

# Solution Strategy
IntegrationScheme                = "Verlet_Velocity"
AutomaticTimestep                = "OFF"
DeltaTimeSafetyFactor            = 1.0
MaxTimeStep                      = 1.0e-7
FinalTime                        = 2.0001
ControlTime                      = 60.0
NeighbourSearchFrequency         = 1

# PostProcess Results
GraphExportFreq                  = 1e-2
VelTrapGraphExportFreq           = 1e-2
OutputTimeStep                   = 1e-2
PostBoundingBox                  = 0
PostDisplacement                 = 1
PostVelocity                     = 1
# DEM only Results
PostTotalForces                  = 0
PostRigidElementForces           = 1
PostRadius                       = 0
PostAngularVelocity              = 0
PostParticleMoment               = 0
PostEulerAngles                  = 1
PostRollingResistanceMoment      = 0
# FEM only Results
PostElasticForces                = 1
PostContactForces                = 1
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
PostStressStrainOption           = 0

#
problem_name="benchmark31"
