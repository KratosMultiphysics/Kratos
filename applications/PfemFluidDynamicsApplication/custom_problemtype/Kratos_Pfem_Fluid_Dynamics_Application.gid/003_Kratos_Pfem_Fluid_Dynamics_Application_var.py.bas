dimension = *GenData(DIMENSION)

#Problem Data
#####################################

ProblemType = "*GenData(Problem_Type)"
NumberofThreads = *GenData(Number_of_threads,INT)
Solution_method = "Newton-Raphson"		
SolverType  = "*GenData(Solver_Type)"
time_step = *GenData(Time_Step)
nsteps = *GenData(Number_of_steps,INT)
EchoLevel = *GenData(Echo_Level)

#Solver Data
#####################################

class SolverSettings:
    echo_level   =  EchoLevel
    solver_type  = "fluid_pfem_mechanical_solver"
    dimension    = *GenData(DIMENSION,INT)
    scheme_type  = "*GenData(Solver_Type)"
*if(strcmp(GenData(Solver_Type),"Dynamic")==0)
    time_integration_method = "*GenData(Time_Integration_Method)"
*else
    time_integration_method = "Implicit"
*endif
*if(strcmp(GenData(Time_Integration_Method),"Explicit")==0)
    explicit_integration_scheme = "*GenData(Explicit_Scheme_Type)"
    max_delta_time = time_step
*endif
*if(strcmp(GenData(Axisymmetric),"True")==0)
    model_type = "Axisymmetric"	    
*else
    model_type = "*GenData(DIMENSION)D"
*endif
*if(strcmp(GenData(DOFS),"ROTATIONS")==0)
    RotationDofs = True
*else
    RotationDofs = False
*endif
*if(strcmp(GenData(DOFS),"U-P")==0)
    PressureDofs = True
    StabilizationFactor = *GenData(Stabilization_Factor)
*else
    PressureDofs = False
*endif
*if(strcmp(GenData(DOFS),"U-wP")==0)
    WaterPressureDofs = True
    StabilizationFactor = *GenData(Stabilization_Factor)
*else
    WaterPressureDofs = False
*endif
*if(strcmp(GenData(Set_Contact_Technology),"1")==0 && strcmp(GenData(Contact_Method),"Lagrange-Multipliers")==0)
    ContactDofs = True
*else
    ContactDofs = False
*endif
    TemperatureDofs = False
*Set cond group_RigidWalls *groups
*if(CondNumEntities > 0)
    RigidWalls = True
*else
    RigidWalls = False
*endif
    ReformDofSetAtEachStep = True
    ComputeReactions = *GenData(Write_Reactions)
    ComputeContactForces = *GenData(Write_Contact_Forces)
    LineSearch = *GenData(LineSearch)
    Implex = *GenData(Implex)
    ComponentWise = *GenData(Component_Wise_Criterion)
    convergence_criterion = "*GenData(Convergence_Criteria)"
*format "%10.5e"
    displacement_relative_tolerance = *GenData(Convergence_Tolerance)
*format "%10.5e"
    displacement_absolute_tolerance = *GenData(Absolute_Tolerance)
*format "%10.5e"		      
    residual_relative_tolerance = *GenData(Convergence_Tolerance)
*format "%10.5e"
    residual_absolute_tolerance = *GenData(Absolute_Tolerance)
*format "%4i"
    max_iteration = *GenData(Max_Iter,INT)

    class linear_solver_config:
	solver_type   = "*GenData(Linear_Solver)"
	scaling       = False
	tolerance     = 1e-7
	max_iteration = *GenData(Linear_Solver_Max_Iteration,INT) #300 
	verbosity     = 0
	is_symmetric  = False  
	
	#Pastix Iterative Solver:
        gmres_krylov_space_dimension = 100
        ilu_level_of_fill            = 3 #5

        #GMRES or CG:
        preconditioner_type          = "None"

        #Deflated CG:
        assume_constant_structure    = True
        max_reduced_size             = 1000

        #AMG: (requires ResidualBasedBlockBuilderAndSolver )
        smoother_type  = "ILU0" #"DAMPED_JACOBI"
        krylov_type    = "GMRES"
	

#Domains Data
#####################################

DomainTypes = []
*set var ndomains(int) = 0
*Set cond group_DeformableBodies *groups
*loop groups *OnlyInCond
*set var ndomains(int)=Operation(ndomains(int)+1)
Domain*ndomains = {"MeshId":*cond(Group_ID), "StructuralType":"*cond(StructuralType)" }
DomainTypes.append(Domain*ndomains)
*end groups


#Meshing Data
#####################################

*set var remesh = 0
*set var refine = 0
*Set cond group_DeformableBodies *groups
*loop groups *OnlyInCond
*if(strcmp(cond(Remesh),"True")==0)
*set var remesh = 1
*else
*set var remesh = 0
*endif
*if(strcmp(cond(Refine),"True")==0)
*set var refine = 1
*else
*set var refine = 0
*endif
*end groups

*if(remesh)
RemeshDomains = True
*else
RemeshDomains = False
*endif
*if(refine)
RefineDomains = True
*else
RefineDomains = False
*endif

*set var ndomains(int) = 0
MeshConditions = []
*loop groups *OnlyInCond
*set var ndomains(int)=Operation(ndomains(int)+1)
string = "*cond(Center_box)" 
BoxCenter*ndomains = [float(s) for s in string.split()]
string = "*cond(Velocity_box)" 
BoxVelocity*ndomains = [float(s) for s in string.split()]
Conditions*ndomains = {"Subdomain":*cond(Group_ID), "StructuralType":"*cond(StructuralType)", "Remesh":*cond(Remesh), "Constrained":*cond(Constrained), "Refine":*cond(Refine), "MeshSmoothing":*cond(MeshSmoothing), "JacobiSmoothing":*cond(JacobiSmoothing), "MeshElement":"*cond(MeshingElement)", "CriticalMeshSize": *cond(Critical_Mesh_Size), "DissipationVariable": "*cond(Dissipation_Variable)", "CriticalDissipation": *cond(Critical_Dissipation), "ErrorVariable": "*cond(Error_Variable)", "CriticalError": *cond(Critical_Error), "TipRadiusRefine": *cond(Tip_Radius_Refine), "CriticalTipRadius": *cond(Critical_Tip_Radius), "RefineOnBoxOnly":*cond(Refine_on_box_only), "BoxCenter": BoxCenter*ndomains, "BoxVelocity": BoxVelocity*ndomains, "BoxRadius": *cond(Radius_box), "RemeshFrequency": *cond(Meshing_Frequency), "LSInterpolation": *cond(LSInterpolation) }
MeshConditions.append(Conditions*ndomains)
*end groups

#set mesh modeler configuration
class mesh_modeler_config:
    echo_level = EchoLevel
    number_domains = *ndomains
    size_scale = 1
    offset_factor = *GenData(Offset_Factor) 
    mesh_conditions = MeshConditions

#Contact Data
#####################################


FindContacts = *GenData(FindContacts)

#set contact modeler configuration
class contact_modeler_config:
    echo_level               = EchoLevel
    contact_condition        = "*GenData(ContactCondition)"
    constrained_contact      = *GenData(Constrained_Contact)
    friction_active          = *GenData(Friction_Active)
    mu_static                = 0.3
    mu_dynamic               = 0.2
    offset_factor            = *GenData(Offset_Factor)
    penalty_parameter        = *GenData(Penalty_Parameter)
    stability_parameter      = *GenData(Stability_Parameter)
    contact_search_frequency = *GenData(Contact_Search_Frequency)


*set var findcontacts = 0
*Set cond group_RigidWalls *groups
*loop groups *OnlyInCond
*if(strcmp(cond(Rigid_Wall_Contact),"True")==0)
*set var findcontacts = 1
*else
*set var findcontacts = 0
*endif
*end groups

*if(findcontacts)
FindRigidWallContacts = True
*else
FindRigidWallContacts = False
*endif

*set var nwalls(int) = 0
WallConditions = []
*loop groups *OnlyInCond
*set var nwalls(int)=Operation(nwalls(int)+1)
*set var nnose(int) = 0
WallNoses*nwalls = [] 
WallPlane = [] 
WallCircle = [] 
*if(strcmp(cond(Wall_Type),"NOSE-WALL")==0)
*for(i=1;i<=cond(Wall_Noses,INT);i=i+7)
*set var nnose(int)=Operation(nnose(int)+1)
tip_center*nnose         = []
tip_center*nnose.append(*cond(Wall_Noses,*i))
tip_center*nnose.append(*cond(Wall_Noses,*Operation(i+1)))
tip_center*nnose.append(*cond(Wall_Noses,*Operation(i+2)))
tip_radius*nnose      = *cond(Wall_Noses,*Operation(i+3))
rake_angle*nnose      = *cond(Wall_Noses,*Operation(i+4))
clearance_angle*nnose = *cond(Wall_Noses,*Operation(i+5))
convexity*nnose       = *cond(Wall_Noses,*Operation(i+6))
Nose*nnose = {"TipCenter": tip_center*nnose, "TipRadius": tip_radius*nnose, "RakeAngle": rake_angle*nnose, "ClearanceAngle": clearance_angle*nnose, "Convexity": convexity*nnose }
WallNoses*nwalls.append(Nose*nnose)
*end
*elseif(strcmp(cond(Wall_Type),"PLANE")==0)
point = []
point.append(*cond(Wall_Plane,1))
point.append(*cond(Wall_Plane,2))
point.append(*cond(Wall_Plane,3))
normal = []
normal.append(*cond(Wall_Plane,4))
normal.append(*cond(Wall_Plane,5))
normal.append(*cond(Wall_Plane,6))
convexity = *cond(Wall_Plane,7)
WallPlane = {"WallPoint": point, "WallNormal": normal, "Convexity": convexity}
*elseif(strcmp(cond(Wall_Type),"CIRCLE")==0)
center = []
center.append(*cond(Wall_Circle,1))
center.append(*cond(Wall_Circle,2))
center.append(*cond(Wall_Circle,3))
radius = *cond(Wall_Circle,4)
convexity = *cond(Wall_Circle,5)
WallCircle = {"Center": center, "Radius": radius, "Convexity": convexity}
*endif

string = "*cond(Linear_Velocity)" 
LinearVelocity*nwalls = [float(s) for s in string.split()]
string = "*cond(Angular_Velocity)" 
AngularVelocity*nwalls = [float(s) for s in string.split()]
string = "*cond(Rotation_Center)" 
RotationCenter*nwalls = [float(s) for s in string.split()]
Conditions*nwalls = {"Subdomain":*cond(Group_ID), "ContactCondition": "*cond(Contact_Condition)", "PenaltyParameter": *cond(Penalty_Parameter), "WallType": "*cond(Wall_Type)", "NumberOfNoses": *nnose, "WallNoses": WallNoses*nwalls, "WallPlane": WallPlane, "WallCircle": WallCircle,"ImposedMovement": *cond(Imposed_Movement), "LinearVelocity": LinearVelocity*nwalls, "AngularVelocity": AngularVelocity*nwalls, "RotationCenter": RotationCenter*nwalls}
WallConditions.append(Conditions*nwalls)
*end groups

#set rigid wall configuration
class rigid_wall_config:
    echo_level         = EchoLevel
    rigid_wall         = FindRigidWallContacts
    size_scale         = 1
    number_of_walls    = *nwalls
    wall_conditions    = WallConditions


#Constraints Data
#####################################

Incremental_Load = "*GenData(Incremental_Load)"
Incremental_Displacement = "*GenData(Incremental_Displacement)"

#PostProcess Data
#####################################

nodal_results = []
nodal_results.append("DISPLACEMENT")
#nodal_results.append("DISPLACEMENT_DT")
*if(strcmp(GenData(Solver_Type),"DynamicSolver")==0)
nodal_results.append("VELOCITY")
nodal_results.append("ACCELERATION")
*endif
*if(strcmp(GenData(Write_Reactions),"True")==0)
nodal_results.append("REACTION")
*endif
*if(strcmp(GenData(Write_Contact_Forces),"True")==0)
nodal_results.append("NORMAL")
nodal_results.append("CONTACT_FORCE")
#nodal_results.append("CONTACT_STRESS")
*endif
*if(strcmp(GenData(DOFS),"U-P")==0)
nodal_results.append("PRESSURE")
*endif
*if(strcmp(GenData(DOFS),"U-wP")==0)
nodal_results.append("WATER_PRESSURE")
*if(strcmp(GenData(Write_Contact_Forces),"True")==0)
#nodal_results.append("EFFECTIVE_CONTACT_FORCE")
#nodal_results.append("EFFECTIVE_CONTACT_STRESS")
*endif
*endif

gauss_points_results = []
*if(strcmp(GenData(Problem_Type),"mechanical")==0)
gauss_points_results.append("GREEN_LAGRANGE_STRAIN_TENSOR")
gauss_points_results.append("CAUCHY_STRESS_TENSOR")
gauss_points_results.append("PLASTIC_STRAIN")
#gauss_points_results.append("INCR_SHEAR_PLASTIC")
#gauss_points_results.append("PRECONSOLIDATION")
#gauss_points_results.append("VOLUMETRIC_PLASTIC")
#gauss_points_results.append("STRESS_INV_P")
#gauss_points_results.append("STRESS_INV_J2")
#gauss_points_results.append("STRESS_INV_THETA")
*endif
*if(strcmp(GenData(DOFS),"U-wP")==0)
#gauss_points_results.append("TOTAL_CAUCHY_STRESS")
#gauss_points_results.append("DARCY_FLOW")
*endif

# GiD output configuration
class GidOutputConfiguration:
    GiDPostMode = "*GenData(File_Format)"
*if(strcmp(GenData(Write_Mesh),"Deformed")==0)
    GiDWriteMeshFlag = True
*else
    GiDWriteMeshFlag = False
*endif
    GiDWriteConditionsFlag = *GenData(Write_Conditions)
    GiDWriteParticlesFlag = *GenData(Write_Particles)
    GiDMultiFileFlag = "Multiples"

GiDWriteFrequency = *GenData(Write_Frequency)
WriteResults = "*GenData(Write_Results)"

# graph_options
PlotGraphs = *GenData(Plot_Graphs)
PlotFrequency = *GenData(Plot_Frequency)

# list options
PrintLists = "*GenData(Print_List_Files)"
file_list = []
*for(i=1;i<=GenData(List_Files,INT);i=i+1)
file_list.append(*GenData(List_Files,*i))
*end

# restart options
SaveRestart = *GenData(Print_Restart)
#RestartFrequency = *GenData(Restart_Frequency)
RestartFrequency=GiDWriteFrequency
LoadRestart = *GenData(Load_Restart)
Restart_Step = *GenData(Load_Step)

### ALL THE WEIGHT THING
TryToSetTheWeight = *GenData(Set_initial_state)
TryToSetConstantWeight = *GenData(Constant_weight)
*if(strcmp(GenData(Set_initial_state),"True")==0)
*if(strcmp(GenData(Constant_weight),"True")==0)
TryToSetWeightVertical = *GenData(SY)
TryToSetWeightHorizontal = *GenData(SX)
*endif
*endif
CPT_PostProcess = *GenData(CPT_PostProcess)


# Declare Python Variables
