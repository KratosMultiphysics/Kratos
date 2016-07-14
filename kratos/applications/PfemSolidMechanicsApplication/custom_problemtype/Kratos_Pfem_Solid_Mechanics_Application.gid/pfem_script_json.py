from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#Activate it to import in the gdb path:
#import sys
#sys.path.append('/home/jmaria/kratos')
#x = input("stopped to allow debug: set breakpoints and press enter to continue");

#### TIME MONITORING START ####

# Time control starts
import time as timer
print(timer.ctime())
# Measure process time
t0p = timer.clock()
# Measure wall time
t0w = timer.time()

def StartTimeMeasuring():
    # Measure process time
    time_ip = timer.clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    # Measure process time
    time_fp = timer.clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[KSM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

#### TIME MONITORING END ####

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication     as KratosSolid
import KratosMultiphysics.PfemBaseApplication           as KratosPfemBase
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid

# Import the python utilities:
import restart_utility                 as restart_utils
import pfem_conditions_python_utility  as condition_utils
import list_files_python_utility       as files_utils
import mesh_modeler_python_utility     as modeler_utils
import rigid_wall_python_utility       as wall_utils
import graph_plot_python_utility       as plot_utils
import time_operation_utility          as operation_utils
import solving_info_utility            as solving_info_utils

######################################################################################
######################################################################################
######################################################################################

#### PARSING THE PARAMETERS ####

#import define_output
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters(parameter_file.read())

#set echo level
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()

#set domain size
domain_size = ProjectParameters["problem_data"]["domain_size"].GetInt()

#### Model_part settings start ####

#defining the model_part
model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : model_part}

#construct the solver (main setting methods are located in the solver_module)
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

# Add variables (always before importing the model part) (it must be integrated in the ImportModelPart)
# If we integrate it in the model part we cannot use combined solvers
solver.AddVariables()

# Read model_part (note: the buffer_size is set here) (restart can be read here)
solver.ImportModelPart()

# Add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
# If we integrate it in the model part we cannot use combined solvers
solver.AddDofs()

# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
## Get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name: main_model_part.GetSubModelPart(part_name)})

#print model_part and properties
if(echo_level>1):
    print("")
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

#### Model_part settings end ####



#### PARSING CLASSICAL PARAMETERS ####


# Import the general variables read from the GiD
import ProjectParameters as general_variables

# defining solver settings
SolverSettings = general_variables.SolverSettings


problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

#TODO: They must be processes:

#### print graph files start ####

plot_active = general_variables.PlotGraphs
graph_plot = plot_utils.GraphPlotUtility(model_part, problem_path)

#### print graph files end ####

#### imposed walls start ####

# set rigid wall contact if it is active:
# activated instead of classical contact
# set rigid wall configuration
rigid_wall = wall_utils.RigidWallUtility(model_part, domain_size, general_variables.rigid_wall_config)

#### imposed walls end ####


#### mesh modeler settings start ####

#construct meshing domains
meshing_domains = []
domains_list = ProjectParameters["meshing_domains"]
for i in range(0,domains_list.size()):
    item = domains_list[i]
    domain_module = __import__(item["python_file_name"].GetString())
    domain = domain_module.CreateMeshingDomain(model_part,item)
    meshing_domains.append(domain)

remesh_domains = general_variables.RemeshDomains
contact_search = general_variables.FindContacts
rigid_wall_contact_search = general_variables.FindRigidWallContacts

modeler = modeler_utils.ModelerUtility(model_part, domain_size, remesh_domains, contact_search, rigid_wall_contact_search)

#### mesh modeler settings end ####


#### manage restart files start ####

# set the restart of the problem
load_restart     = general_variables.LoadRestart
save_restart     = general_variables.SaveRestart

restart_step     = general_variables.Restart_Step
problem_restart  = restart_utils.RestartUtility(model_part, problem_path, problem_name)

#### manage restart files end ####


#set buffer size
buffer_size = 3;

#define problem variables: ( processes add variables )
solver_constructor.AddVariables( model_part, SolverSettings)

# meshing processes
model_part.AddNodalSolutionStepVariable(NORMAL);
model_part.AddNodalSolutionStepVariable(OFFSET);
model_part.AddNodalSolutionStepVariable(SHRINK_FACTOR);
model_part.AddNodalSolutionStepVariable(MEAN_ERROR);
model_part.AddNodalSolutionStepVariable(NODAL_H);
model_part.AddNodalSolutionStepVariable(DETERMINANT_F);

# imposed walls process
model_part.AddNodalSolutionStepVariable(RIGID_WALL);
model_part.AddNodalSolutionStepVariable(WALL_TIP_RADIUS);
model_part.AddNodalSolutionStepVariable(WALL_REFERENCE_POINT);



#--- READ MODEL ------#
if(load_restart == False):
  
  print("::[KPFEM Simulation]:: Reading -START- (MDPA FILE) ")

  #remove results, restart, graph and list previous files
  problem_restart.CleanPreviousFiles()
  list_files.RemoveListFiles()

  #reading the model
  model_part_io = ModelPartIO(problem_name)
  model_part_io.ReadModelPart(model_part)
 
  #set the buffer size
  model_part.SetBufferSize(buffer_size)
  #Note: the buffer size should be set once the mesh is read for the first time

  print("::[KPFEM Simulation]:: Reading -END- ")

  model_part.ProcessInfo[LOAD_RESTART] = 0
    
  #set the degrees of freedom
  solver_constructor.AddDofs(model_part, SolverSettings)

  #set the constitutive law
  import constitutive_law_python_utility as constitutive_law_utils

  constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(model_part,domain_size)
  constitutive_law.Initialize()

else:

  print("::[KPFEM Simulation]:: Reading -RESTART- [FILE",restart_step,"]")

  #reading the model from the restart file
  problem_restart.Load(restart_step);

  print("::[KPFEM Simulation]:: Reading -END- ")

  model_part.ProcessInfo[LOAD_RESTART] = 1

  #remove results, restart, graph and list posterior files
  problem_restart.CleanPosteriorFiles(restart_step)
  list_files.ReBuildListFiles()


# --RIGID WALL OPTIONS START--################
# set rigid wall contact if it is active:
# activated instead of classical contact
# set rigid wall configuration
rigid_wall = wall_utils.RigidWallUtility(model_part, domain_size, general_variables.rigid_wall_config)

# --RIGID WALL OPTIONS END--##################


# --BUILD MESH MODELER START--####################

# build mesh modeler
for domain in meshing_domains:
    domain.SetImposedWalls(rigid_wall)
    domain.Initialize()

modeler.BuildMeshModelers(meshing_domains)

# set rigid walls
#if(rigid_wall_contact_search):
#    modeler.SetRigidWall(rigid_wall)

# --BUILD MESH MODELER END--####################

# --CONTACT SEARCH START--####################

# build mesh modeler
modeler.BuildContactModeler(general_variables.contact_modeler_config)

# --CONTACT SEARCH END--######################


#--- MODELER INITIALIZATION---#

#set mesh searches and modeler
modeler.InitializeDomains( load_restart ); ## due to the skin conditions at reloading

# mesh size nodal h search
if(load_restart == False):
    modeler.SearchNodalH();

#--- PRINT CONTROL ---#

if(echo_level>=1):
    print("")
    print(model_part)
    print(model_part.Properties[1])

# --INITIALIZE--###########################
#

# set delta time in process info
model_part.ProcessInfo[DELTA_TIME] = general_variables.time_step

# solver initialize
main_step_solver.Initialize()
main_step_solver.SetRestart(load_restart) #calls strategy initialize if no restart

# initial contact search
modeler.InitialContactSearch()

#define time steps and loop range of steps
time_step = model_part.ProcessInfo[DELTA_TIME]

if(load_restart):  

  buffer_size  = 0

else:

  model_part.ProcessInfo[TIME]                = 0
  model_part.ProcessInfo[TIME_STEPS]          = 0
  model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step

  conditions.Initialize(time_step);

#initialize step operations
starting_step  = model_part.ProcessInfo[TIME_STEPS]
starting_time  = model_part.ProcessInfo[TIME]
ending_step    = general_variables.nsteps
ending_time    = general_variables.nsteps * time_step


output_print = operation_utils.TimeOperationUtility()
gid_time_frequency = general_variables.GiDWriteFrequency
output_print.InitializeTime(starting_time, ending_time, time_step, gid_time_frequency)

restart_print = operation_utils.TimeOperationUtility()
restart_time_frequency = general_variables.RestartFrequency
restart_print.InitializeTime(starting_time, ending_time, time_step, restart_time_frequency)

mesh_generation = operation_utils.TimeOperationUtility()
mesh_generation_frequency = modeler.GetRemeshFrequency()
mesh_generation.InitializeTime(starting_time, ending_time, time_step, mesh_generation_frequency)

contact_search = operation_utils.TimeOperationUtility()
contact_search_frequency = general_variables.contact_modeler_config.contact_search_frequency
contact_search.InitializeTime(starting_time, ending_time, time_step, contact_search_frequency)

rigid_wall_contact_search = operation_utils.TimeOperationUtility()
rigid_wall_contact_search_frequency = 0
rigid_wall_contact_search.InitializeTime(starting_time, ending_time, time_step, rigid_wall_contact_search_frequency)

#initialize graph plot variables for time integration
if( plot_active == True):
  mesh_id     = 0 #general_variables.PlotMeshId
  x_variable  = "DISPLACEMENT"
  y_variable  = "CONTACT_FORCE"
  graph_plot.Initialize(x_variable,y_variable,mesh_id)

graph_write= operation_utils.TimeOperationUtility()
graph_write_frequency = general_variables.PlotFrequency
graph_write.InitializeTime(starting_time, ending_time, time_step, graph_write_frequency)

solving_print = operation_utils.TimeOperationUtility()
solving_time_frequency = gid_time_frequency
solving_print.InitializeTime(starting_time, ending_time, time_step, solving_time_frequency)



#### Output settings start ####

# Initialize GiD  I/O (gid outputs, file_lists)
from gid_output_process import GiDOutputProcess
output_settings = ProjectParameters["output_configuration"]
gid_output = GiDOutputProcess(computing_model_part,
                              problem_name,
                              output_settings)

gid_output.ExecuteInitialize()

#### Output settings end ####


# --TIME INTEGRATION--#######################
#
  
#writing a single file
gid_print.initialize_results(model_part)

#initialize time integration variables
current_time = starting_time
current_step = starting_step

# filling the buffer
for step in range(0,buffer_size):

  model_part.CloneTimeStep(current_time)
  model_part.ProcessInfo[DELTA_TIME] = time_step
  model_part.ProcessInfo[TIME_STEPS] = step-buffer_size

# writing a initial state results file
current_id = 0
if(load_restart == False):
    if (general_variables.TryToSetTheWeight):
        if (general_variables.TryToSetConstantWeight):
            conditions.SetConstantWeight( general_variables.TryToSetWeightVertical, general_variables.TryToSetWeightHorizontal);
        else:
            conditions.SetWeight();
    gid_print.write_results(model_part, general_variables.nodal_results, general_variables.gauss_points_results, current_time, current_step, current_id)
    list_files.PrintListFiles(current_id);

# set solver info starting parameters
solving_info = solving_info_utils.SolvingInfoUtility(model_part, SolverSettings)

print(" ")
print("::[KPFEM Simulation]:: Analysis -START- ")

# solving the problem
while(current_time < ending_time):

  # store previous time step
  model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step
  # set new time step ( it can change when solve is called )
  time_step = model_part.ProcessInfo[DELTA_TIME]

  current_time = current_time + time_step
  current_step = current_step + 1

  model_part.CloneTimeStep(current_time)
  model_part.ProcessInfo[TIME] = current_time

  # print process information:
  print_info = solving_print.perform_time_operation(current_time)
  if(print_info):
      solving_info.print_step_info(current_time,current_step)
  
  # processes to be executed at the begining of the solution step
  execute_rigid_wall_contact_search = rigid_wall_contact_search.perform_time_operation(current_time)
  if(execute_rigid_wall_contact_search):
      rigid_wall.ExecuteContactSearch()

  #solving the solid problem 
  clock_time = StartTimeMeasuring();

  #solve time step non-linear system
  main_step_solver.Solve()

  StopTimeMeasuring(clock_time,"Solving", False);

  #processes to be executed at the end of the solution step
  rigid_wall.UpdatePosition()

  #plot graphs
  if(plot_active):
    graph_plot.SetStepResult()

  #incremental load
  conditions.SetIncrementalLoad(current_step, time_step);
  ##conditions.CorrectBoundaryConditions(current_step, time_step); ## function to remove load conditions from the contact...
      
  #print the results at the end of the step
  if(general_variables.WriteResults == "PreMeshing"):
    execute_write = output_print.perform_time_operation(current_time)
    if(execute_write):
      clock_time=StartTimeMeasuring();
      current_id = output_print.operation_id()
      #print gid output file
      gid_print.write_results(model_part,general_variables.nodal_results,general_variables.gauss_points_results,current_time,current_step,current_id)
      #print on list files
      list_files.PrintListFiles(current_id);
      solving_info.set_print_info(execute_write, current_id)
      StopTimeMeasuring(clock_time,"Writing Results", False);

  # remesh domains
  execute_meshing = mesh_generation.perform_time_operation(current_time)
  if(execute_meshing):
    modeler.RemeshDomains();
    solving_info.set_meshing_info(execute_meshing,modeler.GetMeshingStep())

  # contact search
  execute_contact_search = contact_search.perform_time_operation(current_time)
  if(execute_contact_search or execute_meshing):
    modeler.ContactSearch();


  # print the results at the end of the step
  if(general_variables.WriteResults == "PostMeshing"):
    execute_write = output_print.perform_time_operation(current_time)
    if(execute_write):
      clock_time=StartTimeMeasuring();
      current_id = output_print.operation_id()
      #print gid output file
      gid_print.write_results(model_part,general_variables.nodal_results,general_variables.gauss_points_results,current_time,current_step,current_id)
      #print on list files
      list_files.PrintListFiles(current_id);
      solving_info.set_print_info(execute_write, current_id)
      StopTimeMeasuring(clock_time,"Writing Results", False);


  # plot graphs
  if(plot_active):
    execute_plot = graph_write.perform_time_operation(current_time)
    if(execute_plot):
      current_id = output_print.operation_id()
      graph_plot.Plot(current_id)

  # print restart file
  if(save_restart):
    execute_save = restart_print.perform_time_operation(current_time)
    if(execute_save):
      clock_time=StartTimeMeasuring();
      current_id = output_print.operation_id()
      problem_restart.Save(current_time,current_step,current_id)
      solving_info.set_restart_info(execute_save,current_id)
      StopTimeMeasuring(clock_time,"Writing Restart", False)


  solving_info.update_solving_info()
  if(print_info):
      solving_info.print_solving_info()

  conditions.RestartImposedDisp()

# --FINALIZE--############################
#

# writing a single file
gid_print.finalize_results()

print("::[KPFEM Simulation]:: Analysis -END- ")
print(" ")

# --END--###############################
#

# check solving information for any problem
solving_info.info_check()

# measure process time
tfp = clock()
# measure wall time
tfw = time()

print("::[KPFEM Simulation]:: [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")

print(ctime())

# to create a benchmark: add standard benchmark files and decomment next two lines 
# rename the file to: run_test.py
#from run_test_benchmark_results import *
#WriteBenchmarkResults(model_part)
