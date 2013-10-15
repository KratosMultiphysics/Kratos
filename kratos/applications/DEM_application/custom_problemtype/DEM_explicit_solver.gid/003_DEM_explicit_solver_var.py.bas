# Model Type

ContinuumOption                  = "*GenData(Continuum)"
RotationOption                   = "*GenData(Rotation)"
HomogeneousMaterialOption        = "*GenData(Homogeneous_Material)"
ElementType                      = "*GenData(DEM_Element_Type)"
Dempack                          = "*GenData(Dempack)"

# Meshing Settings

CleanIndentationsOption          = "*GenData(Clean_Initial_Indentations)"

# General Settings

*format "%10.5e"
FinalTime                        = *GenData(Calculation_Time)
*format "%10.5e"
GravityX                         = *GenData(Gravity_x)
*format "%10.5e"
GravityY                         = *GenData(Gravity_y)
*format "%10.5e"
GravityZ                         = *GenData(Gravity_z)
BoundingBoxOption                = "*GenData(Bounding_Box_option)"
*format "%10.5e"
BoundingBoxEnlargementFactor     = *GenData(Bounding_Box_Enlargement_Factor)
AutomaticBoundingBoxOption       = "*GenData(Automatic_Calculation_Of_Bounding_Box)"
*format "%10.5e"
BoundingBoxMaxX                  = *GenData(Max_X)
*format "%10.5e"
BoundingBoxMaxY                  = *GenData(Max_Y)
*format "%10.5e"
BoundingBoxMaxZ                  = *GenData(Max_Z)
*format "%10.5e"
BoundingBoxMinX                  = *GenData(Min_X)
*format "%10.5e"
BoundingBoxMinY                  = *GenData(Min_Y)
*format "%10.5e"
BoundingBoxMinZ                  = *GenData(Min_Z)
Dimension                        = *GenData(Domain_Dimension)
OutputFileType                   = "*GenData(Output_file_type)"
Multifile                        = "*GenData(Multifile)"
PrintNeighbourLists              = "*GenData(Print_Neighbour_Lists)"
ModelDataInfo                    = "*GenData(Model_Data_Info)"

# Special features

VirtualMassOption                = "*GenData(Virtual_Mass_Option)"
*format "%10.5e"
VirtualMassCoefficient           = *GenData(Virtual_Mass_Coefficient)
*format "%10.5e"
MagicFactor                      = *GenData(Magic_Factor)
DeltaOption                      = "*GenData(Set_Initial_Indentation_To_Zero_Force)"
*format "%10.5e"
SearchRadiusExtension            = *GenData(Search_Radius_Extension)
*format "%10.5e"
AmplifiedSearchRadiusExtension   = *GenData(Amplified_Continuum_Search_Radius_Extension)
FixVelocitiesOption              = "*GenData(Fix_Velocities_At_Predetermined_Time)"
*format "%10.5e"
TotalTimePercentageFixVelocities = *GenData(Time_Step_Constrain_DOFs_Percentage)
TrihedronOption                  = "*GenData(Trihedron_On_Each_Ball)"
LimitSurfaceOption               = *GenData(Limit_Surface)
*format "%10.5e"
SurfaceNormalDirX1               = *GenData(Surface_Normal_Dir_X_1)
*format "%10.5e"
SurfaceNormalDirY1               = *GenData(Surface_Normal_Dir_Y_1)
*format "%10.5e"
SurfaceNormalDirZ1               = *GenData(Surface_Normal_Dir_Z_1)
*format "%10.5e"
SurfacePointCoorX1               = *GenData(Surface_Point_Coor_X_1)
*format "%10.5e"
SurfacePointCoorY1               = *GenData(Surface_Point_Coor_Y_1)
*format "%10.5e"
SurfacePointCoorZ1               = *GenData(Surface_Point_Coor_Z_1)
*format "%10.5e"
SurfaceFrictionAngle1            = *GenData(Surface_Friction_Angle_1)
*format "%10.5e"
SurfaceNormalDirX2               = *GenData(Surface_Normal_Dir_X_2)
*format "%10.5e"
SurfaceNormalDirY2               = *GenData(Surface_Normal_Dir_Y_2)
*format "%10.5e"
SurfaceNormalDirZ2               = *GenData(Surface_Normal_Dir_Z_2)
*format "%10.5e"
SurfacePointCoorX2               = *GenData(Surface_Point_Coor_X_2)
*format "%10.5e"
SurfacePointCoorY2               = *GenData(Surface_Point_Coor_Y_2)
*format "%10.5e"
SurfacePointCoorZ2               = *GenData(Surface_Point_Coor_Z_2)
*format "%10.5e"
SurfaceFrictionAngle2            = *GenData(Surface_Friction_Angle_2)
*format "%10.5e"
SurfaceNormalDirX3               = *GenData(Surface_Normal_Dir_X_3)
*format "%10.5e"
SurfaceNormalDirY3               = *GenData(Surface_Normal_Dir_Y_3)
*format "%10.5e"
SurfaceNormalDirZ3               = *GenData(Surface_Normal_Dir_Z_3)
*format "%10.5e"
SurfacePointCoorX3               = *GenData(Surface_Point_Coor_X_3)
*format "%10.5e"
SurfacePointCoorY3               = *GenData(Surface_Point_Coor_Y_3)
*format "%10.5e"
SurfacePointCoorZ3               = *GenData(Surface_Point_Coor_Z_3)
*format "%10.5e"
SurfaceFrictionAngle3            = *GenData(Surface_Friction_Angle_3)
*format "%10.5e"
SurfaceNormalDirX4               = *GenData(Surface_Normal_Dir_X_4)
*format "%10.5e"
SurfaceNormalDirY4               = *GenData(Surface_Normal_Dir_Y_4)
*format "%10.5e"
SurfaceNormalDirZ4               = *GenData(Surface_Normal_Dir_Z_4)
*format "%10.5e"
SurfacePointCoorX4               = *GenData(Surface_Point_Coor_X_4)
*format "%10.5e"
SurfacePointCoorY4               = *GenData(Surface_Point_Coor_Y_4)
*format "%10.5e"
SurfacePointCoorZ4               = *GenData(Surface_Point_Coor_Z_4)
*format "%10.5e"
SurfaceFrictionAngle4            = *GenData(Surface_Friction_Angle_4)
*format "%10.5e"
SurfaceNormalDirX5               = *GenData(Surface_Normal_Dir_X_5)
*format "%10.5e"
SurfaceNormalDirY5               = *GenData(Surface_Normal_Dir_Y_5)
*format "%10.5e"
SurfaceNormalDirZ5               = *GenData(Surface_Normal_Dir_Z_5)
*format "%10.5e"
SurfacePointCoorX5               = *GenData(Surface_Point_Coor_X_5)
*format "%10.5e"
SurfacePointCoorY5               = *GenData(Surface_Point_Coor_Y_5)
*format "%10.5e"
SurfacePointCoorZ5               = *GenData(Surface_Point_Coor_Z_5)
*format "%10.5e"
SurfaceFrictionAngle5            = *GenData(Surface_Friction_Angle_5)
LimitCylinderOption              = *GenData(Limit_Cylinder)
*format "%10.5e"
CylinderAxisX1                   = *GenData(Cylinder_Axis_X_1)
*format "%10.5e"
CylinderAxisY1                   = *GenData(Cylinder_Axis_Y_1)
*format "%10.5e"
CylinderAxisZ1                   = *GenData(Cylinder_Axis_Z_1)
*format "%10.5e"
CylinderInitialBaseCentreX1      = *GenData(Cylinder_Initial_Base_Centre_X_1)
*format "%10.5e"
CylinderInitialBaseCentreY1      = *GenData(Cylinder_Initial_Base_Centre_Y_1)
*format "%10.5e"
CylinderInitialBaseCentreZ1      = *GenData(Cylinder_Initial_Base_Centre_Z_1)
*format "%10.5e"
CylinderRadius1                  = *GenData(Cylinder_Radius_1)
*format "%10.5e"
CylinderVelocity1                = *GenData(Cylinder_Velocity_1)
*format "%10.5e"
CylinderAngularVelocity1         = *GenData(Cylinder_Angular_Velocity_1)
*format "%10.5e"
CylinderFrictionAngle1           = *GenData(Cylinder_Friction_Angle_1)
*format "%10.5e"
CylinderAxisX2                   = *GenData(Cylinder_Axis_X_2)
*format "%10.5e"
CylinderAxisY2                   = *GenData(Cylinder_Axis_Y_2)
*format "%10.5e"
CylinderAxisZ2                   = *GenData(Cylinder_Axis_Z_2)
*format "%10.5e"
CylinderInitialBaseCentreX2      = *GenData(Cylinder_Initial_Base_Centre_X_2)
*format "%10.5e"
CylinderInitialBaseCentreY2      = *GenData(Cylinder_Initial_Base_Centre_Y_2)
*format "%10.5e"
CylinderInitialBaseCentreZ2      = *GenData(Cylinder_Initial_Base_Centre_Z_2)
*format "%10.5e"
CylinderRadius2                  = *GenData(Cylinder_Radius_2)
*format "%10.5e"
CylinderVelocity2                = *GenData(Cylinder_Velocity_2)
*format "%10.5e"
CylinderAngularVelocity2         = *GenData(Cylinder_Angular_Velocity_2)
*format "%10.5e"
CylinderFrictionAngle2           = *GenData(Cylinder_Friction_Angle_2)
*format "%10.5e"
CylinderAxisX3                   = *GenData(Cylinder_Axis_X_3)
*format "%10.5e"
CylinderAxisY3                   = *GenData(Cylinder_Axis_Y_3)
*format "%10.5e"
CylinderAxisZ3                   = *GenData(Cylinder_Axis_Z_3)
*format "%10.5e"
CylinderInitialBaseCentreX3      = *GenData(Cylinder_Initial_Base_Centre_X_3)
*format "%10.5e"
CylinderInitialBaseCentreY3      = *GenData(Cylinder_Initial_Base_Centre_Y_3)
*format "%10.5e"
CylinderInitialBaseCentreZ3      = *GenData(Cylinder_Initial_Base_Centre_Z_3)
*format "%10.5e"
CylinderRadius3                  = *GenData(Cylinder_Radius_3)
*format "%10.5e"
CylinderVelocity3                = *GenData(Cylinder_Velocity_3)
*format "%10.5e"
CylinderAngularVelocity3         = *GenData(Cylinder_Angular_Velocity_3)
*format "%10.5e"
CylinderFrictionAngle3           = *GenData(Cylinder_Friction_Angle_3)
*format "%10.5e"
CylinderAxisX4                   = *GenData(Cylinder_Axis_X_4)
*format "%10.5e"
CylinderAxisY4                   = *GenData(Cylinder_Axis_Y_4)
*format "%10.5e"
CylinderAxisZ4                   = *GenData(Cylinder_Axis_Z_4)
*format "%10.5e"
CylinderInitialBaseCentreX4      = *GenData(Cylinder_Initial_Base_Centre_X_4)
*format "%10.5e"
CylinderInitialBaseCentreY4      = *GenData(Cylinder_Initial_Base_Centre_Y_4)
*format "%10.5e"
CylinderInitialBaseCentreZ4      = *GenData(Cylinder_Initial_Base_Centre_Z_4)
*format "%10.5e"
CylinderRadius4                  = *GenData(Cylinder_Radius_4)
*format "%10.5e"
CylinderVelocity4                = *GenData(Cylinder_Velocity_4)
*format "%10.5e"
CylinderAngularVelocity4         = *GenData(Cylinder_Angular_Velocity_4)
*format "%10.5e"
CylinderFrictionAngle4           = *GenData(Cylinder_Friction_Angle_4)
*format "%10.5e"
CylinderAxisX5                   = *GenData(Cylinder_Axis_X_5)
*format "%10.5e"
CylinderAxisY5                   = *GenData(Cylinder_Axis_Y_5)
*format "%10.5e"
CylinderAxisZ5                   = *GenData(Cylinder_Axis_Z_5)
*format "%10.5e"
CylinderInitialBaseCentreX5      = *GenData(Cylinder_Initial_Base_Centre_X_5)
*format "%10.5e"
CylinderInitialBaseCentreY5      = *GenData(Cylinder_Initial_Base_Centre_Y_5)
*format "%10.5e"
CylinderInitialBaseCentreZ5      = *GenData(Cylinder_Initial_Base_Centre_Z_5)
*format "%10.5e"
CylinderRadius5                  = *GenData(Cylinder_Radius_5)
*format "%10.5e"
CylinderVelocity5                = *GenData(Cylinder_Velocity_5)
*format "%10.5e"
CylinderAngularVelocity5         = *GenData(Cylinder_Angular_Velocity_5)
*format "%10.5e"
CylinderFrictionAngle5           = *GenData(Cylinder_Friction_Angle_5)

# Time Discretization Settings

IntegrationScheme               = "*GenData(Integration_Scheme)"
*format "%10.5e"
TimeStepsPerSearchStep           = *GenData(Time_Steps_Per_Search_Step)
AutoReductionOfTimeStepOption    = "*GenData(Automatic_Time_Step_Reduction_Option)"
*format "%10.5e"
DeltaTimeSafetyFactor            = *GenData(Time_Step_Safety_Factor)
*format "%10.5e"
MaxTimeStep                      = *GenData(Time_Step)
*format "%10.5e"
OutputTimeStep                   = *GenData(Output_Time_Step)
*format "%10.5e"
ControlTime                      = *GenData(Control_Time)

# Material Model

NormalForceCalculationType       = "*GenData(Normal_Force_Calculation)"
NormalDampingType                = "*GenData(Normal_Contact_Damp)"
TangentialDampingType            = "*GenData(Tangential_Contact_Damp)"
FailureCriterionType             = "*GenData(Failure_Criterion)"
*format "%10.5e"
DempackDamping                   = *GenData(Dempack_Damping)
*format "%10.5e"
TauZero                          = *GenData(Tau_Zero)
*format "%10.5e"
SigmaMax                         = *GenData(Sigma_Max)
*format "%10.5e"
SigmaMin                         = *GenData(Sigma_Min)
*format "%10.5e"
InternalFriction                 = *GenData(Internal_Friction)

*format "%10.5e"
C1                               = *GenData(C_1)
*format "%10.5e"
N1                               = *GenData(N_1)
*format "%10.5e"
C2                               = *GenData(C_2)
*format "%10.5e"
N2                               = *GenData(N_2)
*format "%10.5e"
PlasticYoungModulusRatio         = *GenData(Plastic_Young_Modulus_Ratio)
*format "%10.5e"
PlasticYieldStress               = *GenData(Plastic_Yield_Stress)
*format "%10.5e"
DamageDeformationFactor          = *GenData(Damage_Deformation_Factor)
*format "%10.5e"
G1                               = *GenData(Gamma_1)
*format "%10.5e"
G2                               = *GenData(Gamma_2)
*format "%10.5e"
G3                               = *GenData(Gamma_3)
*format "%10.5e"
MaxDef                           = *GenData(Max_Deformation)

RotationalSpringOption           = "*GenData(Rotation_Spring)"
RotaDampingType                  = "*GenData(Rota_Damp_Type)"

# Global Material Parameters

*format "%10.5e"
GeneralDensity                   = *GenData(General_Density)
*format "%10.5e"
GeneralYoungModulus              = *GenData(General_Young_Modulus)
*format "%10.5e"
GeneralPoissonRatio              = *GenData(General_Poisson_Ratio)
*format "%10.5e"
GeneralCohesion                  = *GenData(General_Cohesion)
*format "%10.5e"
GeneralRollingFriction           = *GenData(General_Rolling_Friction)
*format "%10.5e"
GeneralTension                   = *GenData(General_Tension)
*format "%10.5e"
GeneralRotaDampRatio             = *GenData(General_Rota_Damp_Ratio)
*format "%10.5e"
GeneralStaticFrictionCoef        = *GenData(General_Static_Friction_Coef)
*format "%10.5e"
GeneralDynamicFrictionCoef       = *GenData(General_Dynamic_Friction_Coef)
*format "%10.5e"
GeneralRestitutionCoef           = *GenData(General_Restitution_Coef)
*format "%10.5e"
GeneralColour                    = *GenData(General_Colour)
GlobalVariablesOption            = "*GenData(Globally_Specified_Variables)"
*format "%10.5e"
GlobalKn                         = *GenData(Global_KN)
*format "%10.5e"
GlobalKt                         = *GenData(Global_KT)
*format "%10.5e"
GlobalKr                         = *GenData(Global_KR)
*format "%10.5e"
GlobalRn                         = *GenData(Global_RN)
*format "%10.5e"
GlobalRT                         = *GenData(Global_RT)
*format "%10.5e"
GlobalRr                         = *GenData(Global_RR)
*format "%10.5e"
GlobalFrictionAngle              = *GenData(Global_FRI_ANG)

# Continuum Options

StressStrainOperationsOption     = "*GenData(Stress_Strain_Operations)"
ContactMeshOption                = "*GenData(Contact_Mesh)"
ConcreteTestOption               = "*GenData(Concrete_Test)"
RealTimeGraphOption              = "*GenData(Real_Time_Graph)"
TriaxialOption                   = "*GenData(Triaxial_Option)"
*format "%10.5e"
ConfinementPressure              = *GenData(Confinement_Pressure)
*format "%10.5e"
InitialPressureAplicationTime    = *GenData(Initial_Pressure_Aplication_Time)
*format "%10.5e"
TotalTimePercentAsForceAplTime   = *GenData(Total_Time_Percentage_As_Force_Aplication_Time)

#POSTPROCES

PostVelocity                     = "*GenData(VELOCITY)"
PostDisplacement                 = "*GenData(DISPLACEMENT)"
PostRadialDisplacement           = "*GenData(RADIAL_DISPLACEMENT)"
PostRHS                          = "*GenData(RHS)"
PostTotalForces                  = "*GenData(TOTAL_FORCES)"
PostDampForces                   = "*GenData(DAMP_FORCES)"
PostAppliedForces                = "*GenData(APPLIED_FORCES)"
PostRadius                       = "*GenData(RADIUS)"
PostParticleCohesion             = "*GenData(PARTICLE_COHESION)"
PostParticleTension              = "*GenData(PARTICLE_TENSION)"
PostGroupId                      = "*GenData(GROUP_ID)"
PostExportId                     = "*GenData(EXPORT_ID)"
PostExportParticleFailureId      = "*GenData(EXPORT_PARTICLE_FAILURE_ID)"
PostExportSkinSphere             = "*GenData(EXPORT_SKIN_SPHERE)"
PostLocalContactForce            = "*GenData(LOCAL_CONTACT_FORCE)"
PostFailureCriterionState        = "*GenData(FAILURE_CRITERION_STATE)"
PostContactFailure               = "*GenData(CONTACT_FAILURE)"
PostContactTau                   = "*GenData(CONTACT_TAU)"
PostContactSigma                 = "*GenData(CONTACT_SIGMA)"
PostAngularVelocity              = "*GenData(ANGULAR_VELOCITY)"
PostParticleMoment               = "*GenData(PARTICLE_MOMENT)"
PostEulerAngles                  = "*GenData(EULER_ANGLES)"
PostRepresentativeVolume         = "*GenData(REPRESENTATIVE_VOLUME)"
PostMeanContactArea              = "*GenData(MEAN_CONTACT_AREA)"
PostStressTensor                 = "*GenData(STRESS_TENSOR)"

#FROM CND:

*Set cond volume_SET_SKIN_MANUALLY *elems
*Add cond surface_SET_SKIN_MANUALLY *elems
*Add cond INHERIT_SKIN_FROM_SURFACE *elems
*if(CondNumEntities(int))
PredefinedSkinOption             = "ON"
*else
PredefinedSkinOption             = "OFF"
*endif

TotalElementsVolume              = *tcl(DEM::Get_Mass_Elements)

# Declare Python Variables
