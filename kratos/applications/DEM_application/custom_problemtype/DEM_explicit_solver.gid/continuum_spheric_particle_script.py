import time as timer
import os
import sys
import math
from numpy import *

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import datetime

import DEM_explicit_solver_var as Param
import DEM_procedures
proc = DEM_procedures.Procedures(Param)

import pressure_script as Press

#---------------------MODEL PART KRATOS AND GID.IO ------------------------------------------------------------------

# Defining a model part for the solid part

my_timer = Timer();
balls_model_part = ModelPart("SolidPart");  

# Importing the strategy object

import continuum_sphere_strategy as SolverStrategy

SolverStrategy.AddVariables(balls_model_part, Param)

# Reading the model_part: binary or ascii, multifile or single --> only binary and single for mpi.

if (Param.OutputFileType == "Binary"):
    gid_mode = GiDPostMode.GiD_PostBinary
  
else:
    gid_mode = GiDPostMode.GiD_PostAscii
  
if (Param.Multifile == "multiple_files"):
    multifile = MultiFileFlag.MultipleFiles
  
else:
    multifile = MultiFileFlag.SingleFile

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions   = WriteConditionsFlag.WriteConditions

gid_io = GidIO(Param.problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_solid = ModelPartIO(Param.problem_name)
model_part_io_solid.ReadModelPart(balls_model_part)

# Setting up the buffer size: SHOULD BE DONE AFTER READING!!!

balls_model_part.SetBufferSize(2)

# Adding dofs

SolverStrategy.AddDofs(balls_model_part)

# Constructing a creator/destructor object

creator_destructor = ParticleCreatorDestructor()

# Creating a solver object

solver = SolverStrategy.ExplicitStrategy(balls_model_part, creator_destructor, Param) #here, solver variables initialize as default


# Creating necessary directories

main_path        = os.getcwd()
post_path        = str(main_path) + '/' + str(Param.problem_name) + '_Post_Files'
list_path        = str(main_path) + '/' + str(Param.problem_name) + '_Post_Lists'
neigh_list_path  = str(main_path) + '/' + str(Param.problem_name) + '_Neigh_Lists'
data_and_results = str(main_path) + '/' + str(Param.problem_name) + '_Results_and_Data'
graphs_path      = str(main_path) + '/' + str(Param.problem_name) + '_Graphs'   
MPI_results      = str(main_path) + '/' + str(Param.problem_name) + '_MPI_results'  

for directory in [post_path, list_path, neigh_list_path, data_and_results, graphs_path, MPI_results]:

    if not os.path.isdir(directory):    
        os.makedirs(str(directory))

os.chdir(list_path)

multifile       = open(Param.problem_name + '_all' + '.post.lst', 'w')
multifile_5     = open(Param.problem_name + '_5'   + '.post.lst', 'w')
multifile_10    = open(Param.problem_name + '_10'  + '.post.lst', 'w')
multifile_50    = open(Param.problem_name + '_50'  + '.post.lst', 'w')

multifile.write('Multiple\n')
multifile_5.write('Multiple\n')
multifile_10.write('Multiple\n')
multifile_50.write('Multiple\n')

first_print     = True
index_5         = 1
index_10        = 1
index_50        = 1
prev_time       = 0.0
control         = 0.0


os.chdir(main_path)

export_model_part = balls_model_part

if ( (Param.ContinuumOption =="ON") and (Param.ContactMeshOption =="ON") ) :

  contact_model_part = solver.contact_model_part   
  export_model_part = contact_model_part

#-------------------------------------------------------------------------------------------------------------------------

#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALITZATIONS--------------------------------------------------------

Pressure = 0.0  
Pressure = proc.GiDSolverTransfer(balls_model_part, solver, Param)

if (Param.ModelDataInfo == "ON"):
    os.chdir(data_and_results)
    proc.ModelData(balls_model_part, solver)       # calculates the mean number of neighbours the mean radius, etc..
    os.chdir(main_path)

if (Param.PredefinedSkinOption == "ON" ):

   proc.SetPredefinedSkin(balls_model_part)
 
print 'Initializing Problem....'

dt = Param.MaxTimeStep

total_steps_expected = int(Param.FinalTime / dt)

print ('Total number of TIME STEPs expected in the calculation are: ' + str(total_steps_expected) + ' if time step is kept ' + '\n' )

step_to_fix_velocities = 0

if (Param.FixVelocitiesOption == 'ON'):
  step_to_fix_velocities = 0.01*Param.TotalTimePercentageFixVelocities*total_steps_expected
  print("QUAN HI HAGI CRITICAL TIME STEP... SHA DE MODIFICAR EL TOTAL STEP I EL FIXING..")

balls_model_part.ProcessInfo.SetValue(STEP_FIX_VELOCITIES,int(step_to_fix_velocities))

if( (Param.ContinuumOption == "ON")  and ( (Param.GraphOption =="ON") or (Param.ConcreteTestOption =="ON")) ):
  
    (sup_layer_fm, inf_layer_fm, sup_plate_fm, inf_plate_fm) = proc.ListDefinition(balls_model_part,solver)  # defines the lists where we measure forces

    strain=0.0; total_stress = 0.0; first_time_entry = 1
    # for the graph plotting    
    velocity_node_y = 0.0
    height = 0.3
    diameter = 0.15
  

if(Param.ConcreteTestOption =="ON"):
  
  if(Param.PredefinedSkinOption == "ON" ):
    print "ERROR: in Concrete Test Option the Skin is automatically predefined. Switch the Predefined Skin Option OFF"

  (xtop_area,xbot_area,xlat_area,xtopcorner_area,xbotcorner_area) = proc.CylinderSkinDetermination(balls_model_part,solver,Param) # defines the skin and areas
  
  if ( (Param.TriaxialOption == "ON") and (Pressure != 0.0) ):

    #Correction Coefs
    alpha_top = 3.141592*diameter*diameter*0.25/(xtop_area + 0.70710678*xtopcorner_area)
    alpha_bot = 3.141592*diameter*diameter*0.25/(xbot_area + 0.70710678*xbotcorner_area)
    alpha_lat = 3.141592*diameter*height/(xlat_area + 0.70710678*xtopcorner_area + 0.70710678*xbotcorner_area) 

    Press.ApplyPressure(Pressure, balls_model_part, solver, proc.SKIN, proc.BOT, proc.TOP, proc.LAT, proc.XLAT, proc.XBOT, proc.XTOP, proc.XBOTCORNER, proc.XTOPCORNER, alpha_top, alpha_bot, alpha_lat)

solver.Initialize()

# Initialization of physics monitor and of the initial position of the center of mass

#physics_calculator = SphericElementGlobalPhysicsCalculator(balls_model_part)
#properties_list = []

print 'Initialitzation Complete' + '\n'

os.chdir(graphs_path)

if (Param.GraphOption =="ON"):
  graph_export = open("graph_"+str(datetime.datetime.now())+".csv",'w');
  if (Param.PoissonMeasure =="ON"):
    graph_export_poisson = open("poisson_"+str(datetime.datetime.now())+".csv",'w');

os.chdir(main_path)
    

#Adding stress and strain lists
strainlist=[]; strainlist.append(0.0)
stresslist=[]; stresslist.append(0.0)

if(Param.ContinuumOption =="ON" and Param.GraphOption =="ON"):
  
  #measuring height:

  mean_top = PreUtilities(balls_model_part).MeasureTopHeigh(balls_model_part)
  mean_bot = PreUtilities(balls_model_part).MeasureBotHeigh(balls_model_part)

  ini_height = mean_top - mean_bot
  
  height = ini_height

#
  print ('Initial Height of the Model: ' + str(ini_height)+'\n')


dt1 = dt
dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME)

if(dt != dt1):
  print("sha canviat el pas de temps, corregeix lo de fix vel")

step                   = 0 
time                   = 0.0 
time_old_print         = 0.0
initial_pr_time        = timer.clock()
initial_real_time      = timer.time()
    
#-------------------------------------------------------------------------------------------------------------------------------------

#-----------------------SINGLE FILE MESH AND RESULTS INITIALITZATION-------------------------------------------------------------------

os.chdir(post_path)

if (Param.Multifile == "single_file"):

  if (Param.ContactMeshOption == "ON"):
      gid_io.InitializeMesh(0.0)
      gid_io.WriteMesh(contact_model_part.GetMesh());
      gid_io.FinalizeMesh()
      gid_io.InitializeResults(0.0, contact_model_part.GetMesh()); 

  gid_io.InitializeMesh(0.0)
  gid_io.WriteSphereMesh(balls_model_part.GetMesh())
  gid_io.FinalizeMesh()
  gid_io.InitializeResults(0.0, balls_model_part.GetMesh()); 

#------------------------------------------------------------------------------------------
 
###########################################################################################
#                                                                                         #
#                                    MAIN LOOP                                            #
#                                                                                         #
###########################################################################################
os.chdir(main_path)

print ('Main loop starts at instant: ' + str(initial_pr_time) + '\n')


print ('Total number of TIME STEPs expected in the calculation are: ' + str(total_steps_expected) + ' if time step is kept ' + '\n' )

left_nodes = list()
right_nodes = list()

xleft_weight  = 0.0         
xright_weight  = 0.0

left_counter = 0.0
right_counter = 0.0

if(Param.PoissonMeasure == "ON"):
        
      for node in balls_model_part.Nodes:
        
        if (node.GetSolutionStepValue(GROUP_ID)==4):
          
           left_nodes.append(node)
           xleft_weight = +(node.X0 - node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
           left_counter = +node.GetSolutionStepValue(RADIUS)
           
        elif(node.GetSolutionStepValue(GROUP_ID)==8):
          
           right_nodes.append(node)
           xright_weight = +(node.X + node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
           right_counter = +node.GetSolutionStepValue(RADIUS)
           
      width_ini = xright_weight/right_counter - xleft_weight/left_counter
        

while (time < Param.FinalTime):
 
    dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME
    time = time + dt
    balls_model_part.CloneTimeStep(time)
    balls_model_part.ProcessInfo[TIME_STEPS] = step

    #########################_SOLVE_#########################################4
    os.chdir(main_path) 
    
    solver.Solve()
    #########################TIME CONTROL######################################4
   
    incremental_time = (timer.time() - initial_real_time) - prev_time

    if (incremental_time > Param.ControlTime):
        percentage = 100 * (float(step) / total_steps_expected)
      
        print 'Real time calculation: ' + str(timer.time() - initial_real_time)
        print 'Percentage Completed: '  + str(percentage) + ' %'
        print "TIME STEP = "            + str(step) + '\n'
        
        if( Param.ContinuumOption =="ON" and ( step >= step_to_fix_velocities ) and Param.GraphOption =="ON" and Param.MonitoringOption == "ON"):        
            monitoring = PostUtilities().QuasiStaticAdimensionalNumber(balls_model_part,contact_model_part,balls_model_part.ProcessInfo)
            print "The quasi-static-adimensional-number is:  "            + str(monitoring) + '\n'
            print "The measured stiffness is:  "            + str(total_stress/strain/1e6) + "Mpa" + '\n'
            
        sys.stdout.flush()
        
        prev_time = (timer.time() - initial_real_time)
  
    if ((timer.time() - initial_real_time > 60) and first_print == True and step != 0):    
        first_print = False    
        estimated_sim_duration = 60 * (total_steps_expected / step) # seconds
    
        print('The calculation total estimated time is ' + str(estimated_sim_duration) + 'seconds' + '\n')
        print('in minutes:'        + str(estimated_sim_duration / 60) + 'min.' + '\n')
        print('in hours:'        + str(estimated_sim_duration / 3600) + 'hrs.' + '\n')
        print('in days:'        + str(estimated_sim_duration / 86400) + 'days' + '\n') 
        sys.stdout.flush()
        
        if (estimated_sim_duration / 86400 > 2.0):
          print('WARNING!!!:       VERY LASTING CALCULATION' + '\n')
    
    #########################CONCRETE_TEST_STUFF#########################################4
    
    os.chdir(data_and_results)
    
    total_force = 0.0


    if( Param.ContinuumOption =="ON" and ( step >= step_to_fix_velocities ) and Param.GraphOption =="ON"):
     
      if(first_time_entry):
        #measuring height:

        mean_top = PreUtilities(balls_model_part).MeasureTopHeigh(balls_model_part)
        mean_bot = PreUtilities(balls_model_part).MeasureBotHeigh(balls_model_part)

        ini_height2 = mean_top - mean_bot

        print 'Current Height after confinement: ' + str(ini_height2) + '\n'
        print 'Axial strain due to the confinement: ' + str( 100*(ini_height2-ini_height)/ini_height ) + ' %' +'\n'
        height = ini_height2
        
        for node in sup_layer_fm:
          velocity_node_y = node.GetSolutionStepValue(VELOCITY_Y,0) #Applied velocity during the uniaxial compression test
          break

        print 'velocity for the graph: ' + str(velocity_node_y) + '\n'
              
        first_time_entry = 0

      strain += -2*velocity_node_y*dt/height

      strainlist.append(strain)
      for node in sup_layer_fm:

        force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES,0)[0]
        force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES,0)[1]
        force_node_z = node.GetSolutionStepValue(ELASTIC_FORCES,0)[2]
        
        total_force += force_node_y

      total_stress = total_force/(Param.MeasuringSurface*1000000)
      
      if(Param.PoissonMeasure == "ON"):
                  
        xleft_weight  = 0.0         
        xright_weight  = 0.0

        left_counter = 0.0
        right_counter = 0.0

        for node in left_nodes:
          
          xleft_weight = +(node.X - node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
          left_counter = +node.GetSolutionStepValue(RADIUS)
          
        for node in right_nodes:
          
          xright_weight = +(node.X + node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
          right_counter = +node.GetSolutionStepValue(RADIUS)
        
        width_now = xright_weight/right_counter - xleft_weight/left_counter

        measured_poisson =  ((width_now-width_ini)/width_ini)/strain
        #print( (width_now/0.05)/strain )
    os.chdir(list_path)    
    multifile.write(Param.problem_name + '_' + str(time) + '.post.bin\n')   
    os.chdir(main_path)

  #########################___GiD IO____#########################################4

    time_to_print = time - time_old_print

    if (time_to_print >= Param.OutputTimeStep):
        os.chdir(data_and_results)

        #properties_list = proc.MonitorPhysicalProperties(balls_model_part, physics_calculator, properties_list)

        if (index_5 == 5):       
            multifile_5.write(Param.problem_name + '_' + str(time) + '.post.bin\n')
            index_5 = 0
        
        if (index_10 == 10):       
            multifile_10.write(Param.problem_name + '_' + str(time) + '.post.bin\n')
            index_10 = 0
        
        if (index_50 == 50):
            multifile_50.write(Param.problem_name + '_' + str(time) + '.post.bin\n')
            index_50 = 0
      
        index_5  += 1
        index_10 += 1
        index_50 += 1

        if (Param.Multifile == "multiple_files"):
            gid_io.FinalizeResults()

        os.chdir(graphs_path)

        if (Param.PrintNeighbourLists == "ON"): # Printing neighbours id's
            os.chdir(neigh_list_path)
            neighbours_list = open('neigh_list_' + str(time), 'w')
  
            for elem in balls_model_part.Elements:
                ID = (elem.Id)
                Neigh_ID = elem.GetValue(NEIGHBOURS_IDS)
      
            for i in range(len(Neigh_ID)):
                neighbours_list.write(str(ID) + ' ' + str(Neigh_ID[i]) + '\n')
  
            neighbours_list.close()

        os.chdir(post_path)

        if (Param.Multifile == "multiple_files"):
            gid_io.InitializeMesh(time)
            gid_io.WriteSphereMesh(balls_model_part.GetMesh())
            gid_io.FinalizeMesh()
            gid_io.InitializeResults(time, balls_model_part.GetMesh()); 
 
        proc.PrintingVariables(gid_io, export_model_part, time)
        os.chdir(main_path)     
              
        time_old_print = time
    
    if (Param.GraphOption =="ON"):
      graph_export.write(str(strain)+"  "+str(total_stress)+'\n')
      if (Param.PoissonMeasure =="ON"):
        graph_export_poisson.write(str(strain)+"  "+str(measured_poisson)+'\n')
      
         
    step += 1
#-------------------------------------------------------------------------------------------------------------------------------------


#-----------------------FINALITZATION OPERATIONS-------------------------------------------------------------------------------------- 
#proc.PlotPhysicalProperties(properties_list, graphs_path)

if (Param.Multifile == "single_file"):
    gid_io.FinalizeResults()

graph_export.close()

multifile.close()
multifile_5.close()
multifile_10.close()
multifile_50.close()
os.chdir(main_path)

elapsed_pr_time     = timer.clock() - initial_pr_time
elapsed_real_time   = timer.time() - initial_real_time

print 'Calculation ends at instant: '                 + str(timer.time())
print 'Calculation ends at processing time instant: ' + str(timer.clock())
print 'Elapsed processing time: '                     + str(elapsed_pr_time)
print 'Elapsed real time: '                           + str(elapsed_real_time)

print (my_timer)

print "ANALYSIS COMPLETED" 
