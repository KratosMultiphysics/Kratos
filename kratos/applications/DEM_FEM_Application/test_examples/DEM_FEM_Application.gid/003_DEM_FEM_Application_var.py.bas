# Model Type

ContinuumOption                  = "*GenData(Continuum)"
RotationOption                   = "*GenData(Rotation)"
HomogeneousMaterialOption        = "ON"
ElementType                      = "*GenData(DEM_Element_Type)"

# Meshing Settings

CleanIndentationsOption          = "*GenData(Clean_Initial_Indentations)"

# General Settings

*format "%10.5e"
FinalTime                        = *GenData(TotalTimes)
*format "%10.5e"
GravityX                         = *GenData(Gravity_x)
*format "%10.5e"
GravityY                         = *GenData(Gravity_y)
*format "%10.5e"
GravityZ                         = *GenData(Gravity_z)
BoundingBoxOption                = "OFF"
*format "%10.5e"
BoundingBoxEnlargementFactor     = 1.0
AutomaticBoundingBoxOption       = "OFF"
*format "%10.5e"
BoundingBoxMaxX                  = 1.0
*format "%10.5e"
BoundingBoxMaxY                  = 1.0
*format "%10.5e"
BoundingBoxMaxZ                  = 1.0
*format "%10.5e"
BoundingBoxMinX                  = -1.0
*format "%10.5e"
BoundingBoxMinY                  = -1.0
*format "%10.5e"
BoundingBoxMinZ                  = -1.0
Dimension                        = *GenData(Domain_Dimension)
OutputFileType                   = "*GenData(Output_file_type)"
Multifile                        = "*GenData(Multifile)"
PrintNeighbourLists              = "OFF"
ModelDataInfo                    = "OFF"
FEM_Option                       = "*GenData(FEM_Option)"


# Special features
DeltaOption                      = "*GenData(Search_With_Extension)"
VirtualMassOption                = "*GenData(Virtual_Mass_Option)"
*format "%10.5e"
VirtualMassCoefficient           = *GenData(Virtual_Mass_Coefficient)
*format "%10.5e"
MagicFactor                      = 1.0
SearchRadiusExtension            = *GenData(Search_Radius_Extension)
AmplifiedSearchRadiusExtension   = 1.10000e+00
FixVelocitiesOption              = "OFF"
TotalTimePercentageFixVelocities = 1.00000e+00
TrihedronOption                  = "OFF"
LimitSurfaceOption               = 0
SurfaceNormalDirX1               = 0.00000e+00
SurfaceNormalDirY1               = 1.00000e+00
SurfaceNormalDirZ1               = 0.00000e+00
SurfacePointCoorX1               = 0.00000e+00
SurfacePointCoorY1               = 0.00000e+00
SurfacePointCoorZ1               = 0.00000e+00
SurfaceFrictionAngle1            = 4.50000e+01
SurfaceNormalDirX2               = 0.00000e+00
SurfaceNormalDirY2               = 1.00000e+00
SurfaceNormalDirZ2               = 0.00000e+00
SurfacePointCoorX2               = 0.00000e+00
SurfacePointCoorY2               = 0.00000e+00
SurfacePointCoorZ2               = 0.00000e+00
SurfaceFrictionAngle2            = 4.50000e+01
SurfaceNormalDirX3               = 0.00000e+00
SurfaceNormalDirY3               = 1.00000e+00
SurfaceNormalDirZ3               = 0.00000e+00
SurfacePointCoorX3               = 0.00000e+00
SurfacePointCoorY3               = 0.00000e+00
SurfacePointCoorZ3               = 0.00000e+00
SurfaceFrictionAngle3            = 4.50000e+01
SurfaceNormalDirX4               = 0.00000e+00
SurfaceNormalDirY4               = 1.00000e+00
SurfaceNormalDirZ4               = 0.00000e+00
SurfacePointCoorX4               = 0.00000e+00
SurfacePointCoorY4               = 0.00000e+00
SurfacePointCoorZ4               = 0.00000e+00
SurfaceFrictionAngle4            = 4.50000e+01
SurfaceNormalDirX5               = 0.00000e+00
SurfaceNormalDirY5               = 1.00000e+00
SurfaceNormalDirZ5               = 0.00000e+00
SurfacePointCoorX5               = 0.00000e+00
SurfacePointCoorY5               = 0.00000e+00
SurfacePointCoorZ5               = 0.00000e+00
SurfaceFrictionAngle5            = 4.50000e+01
LimitCylinderOption              = 0
CylinderAxisX1                   = 0.00000e+00
CylinderAxisY1                   = 1.00000e+00
CylinderAxisZ1                   = 0.00000e+00
CylinderInitialBaseCentreX1      = 0.00000e+00
CylinderInitialBaseCentreY1      = 0.00000e+00
CylinderInitialBaseCentreZ1      = 0.00000e+00
CylinderRadius1                  = 1.00000e+00
CylinderVelocity1                = 0.00000e+00
CylinderAngularVelocity1         = 0.00000e+00
CylinderFrictionAngle1           = 4.50000e+01
CylinderAxisX2                   = 0.00000e+00
CylinderAxisY2                   = 1.00000e+00
CylinderAxisZ2                   = 0.00000e+00
CylinderInitialBaseCentreX2      = 0.00000e+00
CylinderInitialBaseCentreY2      = 0.00000e+00
CylinderInitialBaseCentreZ2      = 0.00000e+00
CylinderRadius2                  = 1.00000e+00
CylinderVelocity2                = 0.00000e+00
CylinderAngularVelocity2         = 0.00000e+00
CylinderFrictionAngle2           = 4.50000e+01
CylinderAxisX3                   = 0.00000e+00
CylinderAxisY3                   = 1.00000e+00
CylinderAxisZ3                   = 0.00000e+00
CylinderInitialBaseCentreX3      = 0.00000e+00
CylinderInitialBaseCentreY3      = 0.00000e+00
CylinderInitialBaseCentreZ3      = 0.00000e+00
CylinderRadius3                  = 1.00000e+00
CylinderVelocity3                = 0.00000e+00
CylinderAngularVelocity3         = 0.00000e+00
CylinderFrictionAngle3           = 4.50000e+01
CylinderAxisX4                   = 0.00000e+00
CylinderAxisY4                   = 1.00000e+00
CylinderAxisZ4                   = 0.00000e+00
CylinderInitialBaseCentreX4      = 0.00000e+00
CylinderInitialBaseCentreY4      = 0.00000e+00
CylinderInitialBaseCentreZ4      = 0.00000e+00
CylinderRadius4                  = 1.00000e+00
CylinderVelocity4                = 0.00000e+00
CylinderAngularVelocity4         = 0.00000e+00
CylinderFrictionAngle4           = 4.50000e+01
CylinderAxisX5                   = 0.00000e+00
CylinderAxisY5                   = 1.00000e+00
CylinderAxisZ5                   = 0.00000e+00
CylinderInitialBaseCentreX5      = 0.00000e+00
CylinderInitialBaseCentreY5      = 0.00000e+00
CylinderInitialBaseCentreZ5      = 0.00000e+00
CylinderRadius5                  = 1.00000e+00
CylinderVelocity5                = 0.00000e+00
CylinderAngularVelocity5         = 0.00000e+00
CylinderFrictionAngle5           = 4.50000e+01



# Time Discretization Settings

IntegrationScheme               = "forward_euler"
TimeStepsPerSearchStep           = 1
AutoReductionOfTimeStepOption    = "OFF"
DeltaTimeSafetyFactor            = 1.00000e+00
MaxTimeStep                      = *GenData(Time_Step)
OutputTimeStep                   = *GenData(OutputTimeInterval)
ControlTime                      = *GenData(TotalTimes)


#FROM CND:
PredefinedSkinOption             = "OFF"
TotalElementsVolume              = 0

# Material Model


NormalForceCalculationType       = "*GenData(Normal_Force_Calculation)"
NormalDampingType                = "*GenData(Normal_Contact_Damp)"
TangentialDampingType            = "*GenData(Tangential_Contact_Damp)"



FailureCriterionType             = "Uncoupled"
TauZero                          = 5.26000e+00
SigmaMax                         = 3.32000e+01
SigmaMin                         = 3.32000e+00
InternalFriction                 = 3.50000e+01



RotationalSpringOption           = "*GenData(Rotation_Spring)"
RotaDampingType                  = "*GenData(Rota_Damp_Type)"

C1                               = 4.00000e-01
N1                               = 2.00000e+01
C2                               = 1.50000e+00
N2                               = 1.50000e+01
G1                               = 4.00000e-01
G2                               = 2.00000e+01
G3                               = 1.50000e+00
MaxDef                           = 2.00000e-02

# Global Material Parameters

GeneralDensity                   = 1.00000e+02
GeneralYoungModulus              = 2.10000e+07
GeneralPoissonRatio              = 5.00000e-01
GeneralCohesion                  = 4.16000e+06
GeneralRollingFriction           = 0.00000e+00
GeneralTension                   = 2.01000e+06
GeneralRotaDampRatio             = 5.00000e-01
GeneralStaticFrictionCoef        = 0.00000e+00
GeneralDynamicFrictionCoef       = 0.00000e+00
GeneralRestitutionCoef           = 5.00000e-01
GeneralColour                    = 1.00000e+00
GlobalVariablesOption            = "OFF"
GlobalKn                         = 3.00000e+03
GlobalKt                         = 1.00000e+03
GlobalKr                         = 1.00000e+03
GlobalRn                         = 1.00000e+03
GlobalRT                         = 1.00000e+03
GlobalRr                         = 1.00000e+03
GlobalFrictionAngle              = 4.00000e+01

# Continuum Options

StressStrainOperationsOption     = "OFF"
ContactMeshOption                = "OFF"
ConcreteTestOption               = "OFF"
RealTimeGraphOption              = "OFF"
TriaxialOption                   = "OFF"
ConfinementPressure              = 0.00000e+00
InitialPressureAplicationTime    = 0.00000e+00
TotalTimePercentAsForceAplTime   = 1.50000e+01

#POSTPROCES

PostVelocity                     = "1"
PostDisplacement                 = "1"
PostRadialDisplacement           = "0"
PostRHS                          = "0"
PostTotalForces                  = "0"
PostDampForces                   = "0"
PostAppliedForces                = "0"
PostRadius                       = "1"
PostParticleCohesion             = "0"
PostParticleTension              = "0"
PostGroupId                      = "0"
PostExportId                     = "1"
PostExportParticleFailureId      = "0"
PostExportSkinSphere             = "0"
PostLocalContactForceLow         = "0"
PostLocalContactForceHigh        = "0"
PostFailureCriterionState        = "0"
PostContactFailure               = "0"
PostContactTau                   = "0"
PostContactSigma                 = "0"
PostAngularVelocity              = "0"
PostParticleMoment               = "0"
PostEulerAngles                  = "0"
PostRepresentativeVolume         = "0"
PostMeanContactArea              = "0"
PostStressTensor                 = "0"


ComputeMovementOption            =  "*GenData(ComputeMovementOption)"
RotationSpeed                    = *GenData(RotationSpeed)
AxialSpeed                       = *GenData(AxialSpeed)
PropID                           = *GenData(PropID)
GLOBAL_X_VEL                     = *GenData(GLOBAL_X_VEL)
GLOBAL_Y_VEL                     = *GenData(GLOBAL_Y_VEL)
GLOBAL_Z_VEL                     = *GenData(GLOBAL_Z_VEL)
ROTA_ORIGIN_COORD_X              = *GenData(ROTA_ORIGIN_COORD_X)
ROTA_ORIGIN_COORD_Y              = *GenData(ROTA_ORIGIN_COORD_Y)
ROTA_ORIGIN_COORD_Z              = *GenData(ROTA_ORIGIN_COORD_Z)
ROTA_AXIAL_NORMAL_X              = *GenData(ROTA_AXIAL_NORMAL_X)
ROTA_AXIAL_NORMAL_Y              = *GenData(ROTA_AXIAL_NORMAL_Y)
ROTA_AXIAL_NORMAL_Z              = *GenData(ROTA_AXIAL_NORMAL_Z)
BEGIN_TIME                       = *GenData(BEGIN_TIME)
END_TIME                         = *GenData(END_TIME)




ConstitutiveLaw                  = "*GenData(Constitutive_Law)"
FemDampType                      = "*GenData(FemDampType)"
FemDampRatio                     = *GenData(FemDampRatio)


#FROM CND:

# Declare Python Variables
