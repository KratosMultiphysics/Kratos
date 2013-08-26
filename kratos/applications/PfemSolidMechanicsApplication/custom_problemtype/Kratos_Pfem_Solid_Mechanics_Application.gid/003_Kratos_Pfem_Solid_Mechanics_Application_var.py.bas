domain_size = *GenData(DOMAIN_SIZE)

#Problem Data
#####################################
ProblemType = "*GenData(Problem_Type)"
SolverType  = "*GenData(Solver_Type)"
*format "%10.5e"
time_step   = *GenData(Time_Step)
*format "%8i"
nsteps = *GenData(Number_of_steps,INT)
LineSearch = "*GenData(LineSearch)"
*format "%8i"
NumberofThreads = *GenData(Number_of_threads,INT)

#Solver Data
#####################################
LinearSolver = "*GenData(Linear_Solver)"
*format "%10.5e"
Linear_Solver_Tolerance = *GenData(Linear_Solver_Tolerance)
*format "%8i"
Linear_Solver_Max_Iteration = *GenData(Linear_Solver_Max_Iteration,INT)
Convergence_Criteria = "*GenData(Convergence_Criteria)"
*format "%10.5e"
Convergence_Tolerance = *GenData(Convergence_Tolerance)
*format "%10.5e"
Absolute_Tolerance = *GenData(Absolute_Tolerance)
*format "%4i"
Max_Iter = *GenData(Max_Iter,INT)

#Meshing Data
#####################################
RemeshDomains = "*GenData(RemeshDomains)"
MeshConditions = []
*for(i=1;i<=GenData(MeshDomains,INT);i=i+7)
*set var ndomains(int)=Operation(ndomains(int)+1)
Conditions*ndomains = {"Subdomain":*GenData(MeshDomains,*i), "Remesh":*GenData(MeshDomains,*Operation(i+1)), "Constrained":*GenData(MeshDomains,*Operation(i+2)), "Refine":*GenData(MeshDomains,*Operation(i+3)), "MeshSmoothing":*GenData(MeshDomains,*Operation(i+4)), "JacobiSmoothing":*GenData(MeshDomains,*Operation(i+5)), "MeshElement":"*GenData(MeshDomains,*Operation(i+6))"}
MeshConditions.append(Conditions*ndomains)
*end
NumberDomains  = *ndomains
MeshingElement = "*GenData(MeshingElement)"
RefineDomains  = "False"
MeshingFrequency = *GenData(Meshing_Frequency)
CriticalMeshSize    = *GenData(Critical_Mesh_Size)
CriticalDissipation = *GenData(Critical_Dissipation)
CriticalError       = *GenData(Critical_Error)
RefineBoxOnly       = "*GenData(Refine_on_box_only)"
CenterBox = []
CenterBox.append(*GenData(Center_box,1))
CenterBox.append(*GenData(Center_box,2))
CenterBox.append(0)
VelocityBox = [] 
VelocityBox.append(*GenData(Velocity_box,1))
VelocityBox.append(*GenData(Velocity_box,2))
VelocityBox.append(0)
RadiusBox           = *GenData(Radius_box)

#Contact Data
#####################################
FindContacts = "*GenData(FindContacts)"
ContactCondition = "*GenData(ContactCondition)"
offset_factor = *GenData(Offset_Factor)
penalty_factor = *GenData(Penalty-Stability_Factor)
friction_active = "*GenData(Friction_Active)"
penalty_contact = "*GenData(Penalty_Contact)"
contact_search_frequency = *GenData(Contact_Search_Frequency)
constrained_contact = "*GenData(Constrained_Contact)"
RigidWallContact = "*GenData(Rigid_Wall_Contact)"
tip_radius = *GenData(Tip_Radius)
rake_angle = *GenData(Rake_Angle)
clearance_angle = *GenData(Clearance_Angle)
center = []
center.append(*GenData(Reference_Point_definition,1))
center.append(*GenData(Reference_Point_definition,2))
center.append(0)
velocity = []
velocity.append(*GenData(Velocity_definition,1))
velocity.append(*GenData(Velocity_definition,2))
velocity.append(0)
penalty_parameter = *GenData(Penalty_Parameter)
young_modulus = *GenData(Rigid_Body_Elastic_Modulus)

#Constraints Data
#####################################
Incremental_Load = "*GenData(Incremental_Load)"
Incremental_Displacement = "*GenData(Incremental_Displacement)"

#PostProcess Data
#####################################
echo_level = *GenData(Echo_Level)
WriteResults = "*GenData(WriteResults)"
WriteMesh = "*GenData(WriteMesh)"
FileFormat = "*GenData(FileFormat)"
WriteParticles = "*GenData(WriteParticles)"
WriteConditions = "*GenData(WriteConditions)"
WriteFrequency = *GenData(WriteFrequency)
PlotGraphs = "*GenData(Plot_Graphs)"
PlotFrequency = *GenData(PlotFrequency)
PrintLists = "*GenData(Print_List_Files)"
number_of_lists = *GenData(List_Files)
file_list = []
*for(i=1;i<=GenData(List_Files,INT);i=i+1)
file_list.append(*GenData(List_Files,*i))
*end
SaveRestart = "*GenData(Print_Restart)"
Restart_Interval = *GenData(Restart_Interval)
LoadRestart = "*GenData(Load_Restart)"
Restart_Step = *GenData(Load_Step)

# Declare Python Variables

