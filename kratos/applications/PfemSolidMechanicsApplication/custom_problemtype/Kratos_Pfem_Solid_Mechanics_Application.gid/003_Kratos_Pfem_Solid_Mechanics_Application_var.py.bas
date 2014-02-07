domain_size = *GenData(DOMAIN_SIZE)

#Problem Data
#####################################

ProblemType = "*GenData(Problem_Type)"
NumberofThreads = *GenData(Number_of_threads,INT)
Solution_method = "Newton-Raphson"		
SolverType  = "*GenData(Solver_Type)"
time_step = *GenData(Time_Step)
nsteps = *GenData(Number_of_steps,INT)

#Solver Data
#####################################

class SolverSettings:
    solver_type  = "mechanical_solver"
    domain_size  = *GenData(DOMAIN_SIZE,INT)
    scheme_type = "*GenData(Solver_Type)"
*if(strcmp(GenData(Axisymmetric),"True")==0)
    model_type = "Axisymmetric"	    
*else
    model_type = "*GenData(DOMAIN_SIZE)D"
*endif
*if(strcmp(GenData(DOFS),"ROTATIONS")==0)
    RotationDofs = True
*else
    RotationDofs = False
*endif
*if(strcmp(GenData(DOFS),"U-P")==0)
    PressureDofs = True
*else
    PressureDofs = False
*endif
    TemperatureDofs = False
    RigidWalls = *GenData(Rigid_Wall_Contact)
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
	

#Meshing Data
#####################################

RemeshDomains  = *GenData(RemeshDomains)
MeshingElement = "*GenData(MeshingElement)"
RefineDomains  = False
MeshConditions = []
*for(i=1;i<=GenData(MeshDomains,INT);i=i+7)
*set var ndomains(int)=Operation(ndomains(int)+1)
Conditions*ndomains = {"Subdomain":*GenData(MeshDomains,*i), "Remesh":*GenData(MeshDomains,*Operation(i+1)), "Constrained":*GenData(MeshDomains,*Operation(i+2)), "Refine":*GenData(MeshDomains,*Operation(i+3)), "MeshSmoothing":*GenData(MeshDomains,*Operation(i+4)), "JacobiSmoothing":*GenData(MeshDomains,*Operation(i+5)), "MeshElement":"*GenData(MeshDomains,*Operation(i+6))"}
MeshConditions.append(Conditions*ndomains)
*end

#set mesh modeler configuration
class mesh_modeler_config:
    number_domains = *ndomains
    size_scale = 1
    critical_mesh_size = *GenData(Critical_Mesh_Size)
    critical_dissipation = *GenData(Critical_Dissipation)
    reference_error = *GenData(Critical_Error)
    offset_factor = *GenData(Offset_Factor)
    tip_radius_refine = *GenData(Tip_Radius_Refine)
    critical_tip_radius = *GenData(Critical_Tip_Radius)		  
    mesh_conditions = MeshConditions
    box_refinement_only = *GenData(Refine_on_box_only)
    box_center = []
    box_center.append(*GenData(Center_box,1))
    box_center.append(*GenData(Center_box,2))
    box_center.append(0)
    box_velocity = []
    box_velocity.append(*GenData(Velocity_box,1))
    box_velocity.append(*GenData(Velocity_box,2))
    box_velocity.append(0)
    box_radius = *GenData(Radius_box)
    remesh_frequency = *GenData(Meshing_Frequency)

#Contact Data
#####################################

FindContacts = *GenData(FindContacts)
FindRigidWallContacts = *GenData(Rigid_Wall_Contact)

class contact_modeler_config:
    contact_condition        = "*GenData(ContactCondition)"
    constrained_contact      = *GenData(Constrained_Contact)
    friction_active          = *GenData(Friction_Active)
    mu_static                = 0.3
    mu_dynamic               = 0.2
    offset_factor            = *GenData(Offset_Factor)
    penalty_parameter        = *GenData(Penalty_Parameter)
    stability_parameter      = *GenData(Stability_Parameter)
    contact_search_frequency = *GenData(Contact_Search_Frequency)

#set rigid wall configuration
class rigid_wall_config:
    rigid_wall         = *GenData(Rigid_Wall_Contact)
    size_scale         = 1
    contact_condition  = "*GenData(Contact_Condition)"
    penalty_parameter  = *GenData(Rigid_Body_Penalty_Parameter)	
    number_of_walls    = *GenData(Number_of_walls)
    wall_labels        = []
    tip_radius         = []
    rake_angles        = []   	 
    clearance_angles   = []  
    tip_centers        = []
    nose_convexities   = []     
*for(i=1;i<=GenData(Wall_Noses,INT);i=i+8)
    wall_labels.append(*GenData(Wall_Noses,*i))  
    tip_radius.append(*GenData(Wall_Noses,*Operation(i+1)))
    rake_angles.append(*GenData(Wall_Noses,*Operation(i+2)))
    clearance_angles.append(*GenData(Wall_Noses,*Operation(i+3)))
    tip_centers.append([*GenData(Wall_Noses,*Operation(i+4)),*GenData(Wall_Noses,*Operation(i+5)),*GenData(Wall_Noses,*Operation(i+6))])
    nose_convexities.append(*GenData(Wall_Noses,*Operation(i+7)))
*end
    wall_movement_labels  = [] 
    wall_velocity         = []
*for(i=1;i<=GenData(Wall_Movements,INT);i=i+4)
    wall_movement_labels.append(*GenData(Wall_Movements,*i))
    wall_velocity.append([*GenData(Wall_Movements,*Operation(i+1)),*GenData(Wall_Movements,*Operation(i+2)),*GenData(Wall_Movements,*Operation(i+3))])
*end
    wall_rotation_labels  = [] 
    wall_angular_velocity = []
    reference_point       = []
*for(i=1;i<=GenData(Wall_Rotations,INT);i=i+7)
    wall_rotation_labels.append(*GenData(Wall_Rotations,*i))
    wall_angular_velocity.append([*GenData(Wall_Rotations,*Operation(i+1)),*GenData(Wall_Rotations,*Operation(i+2)),*GenData(Wall_Rotations,*Operation(i+3))]) 
    reference_point.append([*GenData(Wall_Rotations,*Operation(i+4)),*GenData(Wall_Rotations,*Operation(i+5)),*GenData(Wall_Rotations,*Operation(i+6))])     
*end

#Constraints Data
#####################################

Incremental_Load = "*GenData(Incremental_Load)"
Incremental_Displacement = "*GenData(Incremental_Displacement)"

#PostProcess Data
#####################################

nodal_results = []
nodal_results.append("DISPLACEMENT")
*if(strcmp(GenData(Write_Reactions),"True")==0)
nodal_results.append("REACTION")
*endif
*if(strcmp(GenData(Write_Contact_Forces),"True")==0)
nodal_results.append("NORMAL")
nodal_results.append("CONTACT_FORCE")
*endif

gauss_points_results = []
gauss_points_results.append("GREEN_LAGRANGE_STRAIN_TENSOR")
gauss_points_results.append("CAUCHY_STRESS_TENSOR")
gauss_points_results.append("PLASTIC_STRAIN")

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
echo_level = *GenData(Echo_Level)

# graph_options
PlotGraphs = "*GenData(Plot_Graphs)"
PlotFrequency = *GenData(Plot_Frequency)

# list options
PrintLists = "*GenData(Print_List_Files)"
file_list = []
*for(i=1;i<=GenData(List_Files,INT);i=i+1)
file_list.append(*GenData(List_Files,*i))
*end

# restart options
SaveRestart = *GenData(Print_Restart)
RestartFrequency = *GenData(Restart_Frequency)
LoadRestart = *GenData(Load_Restart)
Restart_Step = *GenData(Load_Step)

# Declare Python Variables

