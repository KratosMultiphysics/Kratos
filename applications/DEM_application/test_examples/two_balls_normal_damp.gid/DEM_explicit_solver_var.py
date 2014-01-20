from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Model Type

ContinuumOption = "OFF"
RotationOption = "OFF"
HomogeneousMaterialOption = "OFF"
ElementType = "SphericParticle3D"

# Meshing Settings

CleanIndentationsOption = "OFF"

# General Settings

FinalTime = 1.00000e+00
GravityX = 0.00000e+00
GravityY = -9.81000e+00
GravityZ = 0.00000e+00
BoundingBoxOption = "OFF"
BoundingBoxEnlargementFactor = 2.00000e+00
Dimension = 3
OutputFileType = "Binary"
Multifile = "single_file"
PrintNeighbourLists = "OFF"
ModelDataInfo = "OFF"

# Special features

VirtualMassOption = "OFF"
VirtualMassCoefficient = 0.00000e+00
MagicFactor = 1.00000e+00
DeltaOption = "OFF"
SearchRadiusExtension = 1.00000e-02
AmplifiedSearchRadiusExtension = 1.10000e+00
FixVelocitiesOption = "OFF"
TotalTimePercentageFixVelocities = 1.00000e+01
TrihedronOption = "OFF"
LimitSurfaceOption = "OFF"
SurfaceNormalDirX = 0.00000e+00
SurfaceNormalDirY = 1.00000e+00
SurfaceNormalDirZ = 0.00000e+00
SurfacePointCoorX = 0.00000e+00
SurfacePointCoorY = 0.00000e+00
SurfacePointCoorZ = 0.00000e+00
SurfaceFrictionAngle = 45
CylinderRadius = 1.00000e+00

# Time Discretization Settings

IntegrationScheme = "forward_euler"
TimeStepsPerSearchStep = 1.00000e+00
AutoReductionOfTimeStepOption = "ON"
DeltaTimeSafetyFactor = 1.00000e+00
MaxTimeStep = 1.00000e-04
OutputTimeStep = 1.00000e-02
ControlTime = 2.00000e+01

# Material Model

NormalForceCalculationType = "Linear"
NormalDampingType = "ViscDamp"
TangentialDampingType = "NoDamp"
FailureCriterionType = "Uncoupled"
TauZero = 5.26000e+00
SigmaMax = 3.32000e+01
SigmaMin = 3.32000e+00
InternalFriction = 3.50000e+01
NonLinearNormalElasticOption = "OFF"
C1 = 4.00000e-01
N1 = 2.00000e+01
C2 = 1.50000e+00
N2 = 1.50000e+01
RotationalSpringOption = "OFF"
RotaDampingType = "NoDamp"

# Global Material Parameters

GeneralDensity = "1.00000e+02"
GeneralYoungModulus = "2.10000e+07"
GeneralPoissonRatio = "5.00000e-01"
GeneralCohesion = "4.16000e+06"
GeneralRollingFriction = "0.00000e+00"
GeneralTension = "2.01000e+06"
GeneralRotaDampRatio = "5.00000e-01"
GeneralStaticFrictionCoef = "0.00000e+00"
GeneralDynamicFrictionCoef = "0.00000e+00"
GeneralRestitutionCoef = "5.00000e-01"
GeneralColour = "1.00000e+00"
GlobalVariablesOption = "OFF"
GlobalKn = 3.00000e+03
GlobalKt = 1.00000e+03
GlobalKr = 1.00000e+03
GlobalRn = 1.00000e+03
GlobalRT = 1.00000e+03
GlobalRr = 1.00000e+03
GlobalFrictionAngle = 4.00000e+01

# Continuum Options

StressStrainOperationsOption = "OFF"
ContactMeshOption = "OFF"
ConcreteTestOption = "OFF"
RealTimeGraphOption = "OFF"
TriaxialOption = "OFF"
ConfinementPressure = 0.00000e+00
InitialPressureAplicationTime = 0.00000e+00
TotalTimePercentAsForceAplTime = 1.50000e+01

# POSTPROCES

PostVelocity = "1"
PostDisplacement = "1"
PostRadialDisplacement = "0"
PostRHS = "0"
PostTotalForces = "0"
PostDampForces = "0"
PostAppliedForces = "0"
PostRadius = "0"
PostParticleCohesion = "0"
PostParticleTension = "0"
PostGroupId = "0"
PostExportId = "0"
PostExportParticleFailureId = "0"
PostExportSkinSphere = "0"
PostLocalContactForceLow = "0"
PostLocalContactForceHigh = "0"
PostFailureCriterionState = "0"
PostContactFailure = "0"
PostContactTau = "0"
PostContactSigma = "0"
PostAngularVelocity = "0"
PostParticleMoment = "0"
PostEulerAngles = "0"
PostRepresentativeVolume = "0"
PostMeanContactArea = "0"
PostStressTensor = "0"

# FROM CND:

PredefinedSkinOption = "OFF"

TotalElementsVolume = 2.26195e-01

# Declare Python Variables
problem_name = 'two_balls_normal_damp'
problem_path = '/home/cimne/kratos/applications/DEM_application/test_examples/two_balls_normal_damp.gid'
kratos_path = '/home/cimne/kratos'
