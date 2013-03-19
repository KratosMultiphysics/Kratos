import time as timer
import os
import sys
import math

#
#
#

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.mpi import *

from DEM_explicit_solver_var import *
from DEM_procedures import *

#---------------------MODEL PART KRATOS AND GID.IO ------------------------------------------------------------------

#defining a model part for the solid part
my_timer=Timer();
solid_model_part = ModelPart("SolidPart");  

import sphere_strategy as SolverStrategy
SolverStrategy.AddVariables(solid_model_part)

AddMpiVariables(solid_model_part)

## reading the solid part: binary or ascii, multifile or single --> only binary and single for mpi.

#
#
#
gid_mode = GiDPostMode.GiD_PostBinary
#
#
#
#
multifile = MultiFileFlag.SingleFile

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions

gid_io = GidIO(problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_solid = ModelPartIO(problem_name)
model_part_io_solid = PerformInitialPartition(solid_model_part,model_part_io_solid,problem_name)
model_part_io_solid.ReadModelPart(solid_model_part)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
solid_model_part.SetBufferSize(2)

##adding dofs
SolverStrategy.AddDofs(solid_model_part)
print('\n')

#creating a solver object

solver = SolverStrategy.ExplicitStrategy(solid_model_part, domain_size) #here, solver variables initialize as default

#----------------------------------------------------------------------------------------------------------------------

#----------------------PYTHON STUFF:------------------------------------------------------------------------------------

#main_path    = os.getcwd()
#post_path    = str(main_path)+'/'+str(problem_name)+'_Post_Files'
#list_path    = str(main_path)+'/'+str(problem_name)+'_Post_Lists'
#neigh_list_path  = str(main_path)+'/'+str(problem_name)+'_Neigh_Lists'
#data_and_results = str(main_path)+'/'+str(problem_name)+'_Results_and_Data'
#graphs_path  = str(main_path)+'/'+str(problem_name)+'_Graphs'   
#MPI_results    = str(main_path)+'/'+str(problem_name)+'_MPI_results'    

#

#
    
#

#

#
#

#
prev_time = 0.0; control = 0.0; cond = 0

#

export_model_part = solid_model_part

if ( (ContinuumOption =="ON") and (ContactMeshOption =="ON") ) :
  
  contact_model_part = solver.contact_model_part   
  export_model_part = contact_model_part
 
#-------------------------------------------------------------------------------------------------------------------------

#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALITZATION--------------------------------------------------------

Pressure = 0.0  
Pressure = ProcGiDSolverTransfer(solid_model_part,solver)

if(ModelDataInfo =="ON"):
#
  ProcModelData(solid_model_part,solver)       # calculates the mean number of neighbours the mean radius, etc..
#
  
if (predefined_skin_option == "ON" ):
   
   for element in solid_model_part.Elements:
         
      element.SetValue(SKIN_SPHERE,0)   
   
      if (element.GetValue(PREDEFINED_SKIN)>0.0): #PREDEFINED_SKIN is a double
      
         element.SetValue(SKIN_SPHERE,1)

if(mpi.rank == 0):
  print'Initialitzating Problem....'
solver.Initialize()
if(mpi.rank == 0):
  print 'Initialitzation Complete' + '\n'

if(ConcreteTestOption =="ON"):
  (sup_layer_fm, inf_layer_fm, sup_plate_fm, inf_plate_fm) = ProcListDefinition(solid_model_part,solver)  # defines the lists where we measure forces

  (xtop_area,xbot_area,xlat_area,xtopcorner_area,xbotcorner_area) = ProcSkinAndPressure(solid_model_part,solver) # defines the skin and areas

  strain=0.0; total_stress = 0.0; first_time_entry = 1
  # for the graph plotting    
  velocity_node_y = 0.0
  height = 0.3
  diameter = 0.15
  
  if ( (TriaxialOption == "ON") and (Pressure != 0.0) ):
      
    xtop_area_gath            = mpi.allgather(mpi.world, xtop_area) 
    xbot_area_gath            = mpi.allgather(mpi.world, xbot_area) 
    xlat_area_gath            = mpi.allgather(mpi.world, xlat_area) 
    xtopcorner_area_gath      = mpi.allgather(mpi.world, xtopcorner_area) 
    xbotcorner_area_gath      = mpi.allgather(mpi.world, xbotcorner_area) 
    
    xtop_area = reduce(lambda x,y:x+y, xtop_area_gath)
    xbot_area = reduce(lambda x,y:x+y, xbot_area_gath)
    xlat_area = reduce(lambda x,y:x+y, xlat_area_gath)
    xtopcorner_area = reduce(lambda x,y:x+y, xtopcorner_area_gath)
    xbotcorner_area = reduce(lambda x,y:x+y, xbotcorner_area_gath)
    
    #Correction Coefs
    alpha_top = 3.141592*diameter*diameter*0.25/(xtop_area + 0.70710678*xtopcorner_area)
    alpha_bot = 3.141592*diameter*diameter*0.25/(xbot_area + 0.70710678*xbotcorner_area)
    alpha_lat = 3.141592*diameter*height/(xlat_area + 0.70710678*xtopcorner_area + 0.70710678*xbotcorner_area) 
      
    ProcApplyPressure(Pressure,solid_model_part,solver,alpha_top,alpha_bot,alpha_lat)
    
if(mpi.rank == 0):
  graph_export = open("strain_stress_data.csv",'w');

#Adding stress and strain lists
#
#

if(ContinuumOption =="ON" and ConcreteTestOption =="ON"):
  
  (Y_mean_bot,counter_bot) = ProcMeasureBOT(BOT,solver)
  (Y_mean_top,counter_top) = ProcMeasureTOP(TOP,solver)
  
  Y_mean_bot_gath   = mpi.gather(mpi.world, Y_mean_bot, 0) 
  counter_bot_gath  = mpi.gather(mpi.world, counter_bot, 0) 
  Y_mean_top_gath   = mpi.gather(mpi.world, Y_mean_top, 0) 
  counter_top_gath  = mpi.gather(mpi.world, counter_top, 0)
  
  if(mpi.rank == 0):
      Y_mean_bot = reduce(lambda x,y:x+y, Y_mean_bot_gath)
      counter_bot = reduce(lambda x,y:x+y, counter_bot_gath)
      Y_mean_top = reduce(lambda x,y:x+y, Y_mean_top_gath)
      counter_top = reduce(lambda x,y:x+y, counter_top_gath)
  
      ini_height = Y_mean_top/counter_top - Y_mean_bot/counter_bot
      
      height = ini_height
      if(mpi.rank==0): 
        print ('Initial Height of the Model: ' + str(ini_height)+'\n')

dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

if (CriticalTimeOption =="ON"):
  solver.Initial_Critical_Time() 

  if (dt!=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)):
    dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)
    if(mpi.rank==0): 
      print("WARNING: Delta time has been modifyed to the critical one")   

time = 0.0; step = 0; time_old_print = 0.0

initial_pr_time = timer.clock()
initial_real_time = timer.time()

if(mpi.rank==0): 
  print ('SOLVE starts at instant: ' + str(initial_pr_time)+'\n')

total_steps_expected = int(final_time/dt)
if(mpi.rank==0): 
  print ('Total number of TIME STEPs expected in the calculation are: ' + str(total_steps_expected) + ' if time step is kept ' +'\n' )
    
#-------------------------------------------------------------------------------------------------------------------------------------

#-----------------------SINGLE FILE MESH AND RESULTS INITIALITZATION-------------------------------------------------------------------

#

gid_io.ChangeOutputName(problem_name + "_" + str(mpi.rank))

if (ContactMeshOption =="ON"): 
    gid_io.InitializeMesh(0.0)
    gid_io.WriteMesh(contact_model_part.GetMesh());
    gid_io.FinalizeMesh()
    gid_io.InitializeResults(0.0, contact_model_part.GetMesh()); 

gid_io.InitializeMesh(0.0)
gid_io.WriteSphereMesh(solid_model_part.GetMesh())
gid_io.FinalizeMesh()
gid_io.InitializeResults(0.0, solid_model_part.GetMesh()); 

#------------------------------------------------------------------------------------------
 
###########################################################################################
#                                                                                         #
#                                    MAIN LOOP                                            #
#                                                                                         #
###########################################################################################
#
while(time < final_time):
 
    dt = solid_model_part.ProcessInfo.GetValue(DELTA_TIME) #possible modifications of DELTA_TIME
    time = time + dt
    solid_model_part.CloneTimeStep(time)

    solid_model_part.ProcessInfo[TIME_STEPS] = step

    #########################_SOLVE_#########################################4
#
    solver.Solve()
    #########################TIME CONTROL######################################4
   
    incremental_time = (timer.time()-initial_real_time)- prev_time

    if (  (incremental_time > control_time) and (mpi.rank == 0) ): 
      
      percentage = 100.0*(float(step)/total_steps_expected)
      print 'Real time calculation: ' + str(timer.time()-initial_real_time) 
      print 'Percentage Completed: ' +str(percentage) + ' %' 
      print "TIME STEP = " + str(step) + '\n'
        
      prev_time = (timer.time()-initial_real_time)
  
    if ( (timer.time()-initial_real_time > 60.0) and cond==0 and mpi.rank==0 ):
    
      cond=1
    
      estimation_time=60.0*(total_steps_expected/step) #seconds
    
      print('the total calculation estimated time is '+str(estimation_time)+'seconds.'+'\n')
      print('in minutes :'+str(estimation_time/60)+'min.'+'\n')
      print('in hours :'+str((estimation_time/60)/60)+'hrs.'+'\n')
      print('in days :'+str(((estimation_time/60)/60)/24)+'days.'+'\n') 
    
      if (((estimation_time/60)/60)/24 > 2.0):
        print('WARNING!!!:       VERY LASTING CALCULATION'+'\n')
    
    #########################CONCRETE_TEST_STUFF#########################################4
 
    if( (ConcreteTestOption =="ON") and (step==2) ):
      
      #Cross section Area Control
      if(mpi.rank==0): 
        Num_Cross_Sect = solid_model_part.ProcessInfo.GetValue(AREA_VERTICAL_TAPA)
        Exact_Cross_Sect = 3.141592*0.15*0.15*0.25
        
        print '\n' + '----------------------CONCRETE TEST CONTROLS----------------------' + '\n'
        print 'Total Horitzontal Numerical Cross Section on Force Measurement: ' + str(Num_Cross_Sect)
        print 'Total Horitzontal Real Cross Section on Force Measurement was: ' + str(Exact_Cross_Sect)
        print 'Relative Error: ' + str (100*(abs(Num_Cross_Sect-Exact_Cross_Sect)/Exact_Cross_Sect)) + ' %'
        #print( solid_model_part.ProcessInfo.GetValue(AREA_VERTICAL_CENTRE) )
      
      total_volume = 0.0;  h   = 0.3;    d   = 0.15
      
      for element in solid_model_part.Elements:
  
        node = element.GetNode(0)
        volume_equ = node.GetSolutionStepValue(REPRESENTATIVE_VOLUME,0) 
        total_volume += volume_equ
     
      total_volume_gath = mpi.gather(mpi.world, total_volume, 0)
      
      if(mpi.rank==0):       
        total_volume = reduce(lambda x,y:x+y, total_volume_gath)
        real_volume = 3.141592*d*d*0.25*h     
        print 'Total Numerical Volume: ' + str(total_volume)
        print 'Total Real Volume: ' + str(real_volume)
        print 'Error: ' + str(100*abs(total_volume-real_volume)/real_volume) +'%'+'\n'
        print '------------------------------------------------------------------' + '\n'
    
#
    
    total_force=0.0
    force_node= 0.0
    
    #For a uniaxial compression test with a cylinder of 15 cm diameter and 30 cm height

    if( FixVelocities == 'OFF'):
      TimePercentageFixVelocities = 0.0
      
    if( ContinuumOption =="ON" and ( time >= 0.01*TimePercentageFixVelocities*final_time) and ConcreteTestOption =="ON"):
     
      if(first_time_entry):
        (Y_mean_bot,counter_bot) = ProcMeasureBOT(BOT,solver)
        (Y_mean_top,counter_top) = ProcMeasureTOP(TOP,solver)
       
        Y_mean_bot_gath   = mpi.gather(mpi.world, Y_mean_bot, 0) 
        counter_bot_gath   = mpi.gather(mpi.world, counter_bot, 0) 
        Y_mean_top_gath   = mpi.gather(mpi.world, Y_mean_top, 0) 
        counter_top_gath   = mpi.gather(mpi.world, counter_top, 0)
  
        if(mpi.rank == 0):
          Y_mean_bot = reduce(lambda x,y:x+y, Y_mean_bot_gath)
          counter_bot = reduce(lambda x,y:x+y, counter_bot_gath)
          Y_mean_top = reduce(lambda x,y:x+y, Y_mean_top_gath)
          counter_top = reduce(lambda x,y:x+y, counter_top_gath)
        
          ini_height2 = Y_mean_top/counter_top - Y_mean_bot/counter_bot
         
          print 'Current Height after confinement: ' + str(ini_height2) + '\n'
          print 'Axial strain due to the confinement: ' + str( 100*(ini_height2-ini_height)/ini_height ) + ' %' +'\n'
          height = ini_height2
          
        for node in sup_layer_fm:
          velocity_node_y = node.GetSolutionStepValue(VELOCITY_Y,0) #Applied velocity during the uniaxial compression test
          break
        
        velocity_gath   = mpi.gather(mpi.world, velocity_node_y, 0) 
        
        if(mpi.rank == 0):
          for vel in velocity_gath:
            if (vel != 0.0):
              velocity_node_y = vel              #only if all are the same
              print 'velocity for the graph: ' + str(velocity_node_y) + '\n'
              break
              
        first_time_entry = 0

      strain += -2*velocity_node_y*dt/height

      #

      for node in sup_layer_fm:
      
        force_node = node.GetSolutionStepValue(RHS,0)
        force_node_x = node.GetSolutionStepValue(RHS,0)[0]
        force_node_y = node.GetSolutionStepValue(RHS,0)[1]
        force_node_z = node.GetSolutionStepValue(RHS,0)[2]
        
        total_force += force_node_y
            
      total_force_gath   = mpi.gather(mpi.world, total_force, 0) 
      if(mpi.rank == 0):
        total_force = reduce(lambda x,y:x+y,total_force_gath)
        total_stress = total_force/(math.pi*75*75) + (1e-6)*Pressure #Stress in MPa
    #

    #
    
    #
    
    #os.chdir(main_path)

  #########################___GiD IO____#########################################4

    time_to_print = time - time_old_print

    if(time_to_print >= output_dt):

      #
        
      #
        
      #
      #
        
      #
        
      #
      #
        
      #
      #
      #
      
      #

      #
      #

      #

      #

      #
        
      #
      #
      #
      #
      #
      #
      #

      if(PrintNeighbourLists == "ON"): #printing neighbours id's
  
        #os.chdir(neigh_list_path)
        neighbours_list = open('neigh_list_'+ str(time),'w')
  
        for elem in solid_model_part.Elements:
          ID=(elem.Id)
          Neigh_ID = elem.GetValue(NEIGHBOURS_IDS)
      
          for i in range(len(Neigh_ID)):

            neighbours_list.write(str(ID)+' '+str(Neigh_ID[i])+'\n')
  
        neighbours_list.close()

      #os.chdir(post_path)

      #

      #
      #
      #
      #
      
      #
      #
      #
      #
      #
 
      ProcPrintingVariables(gid_io,export_model_part,time)  

      #os.chdir(main_path)     
              
      time_old_print = time

    if(mpi.rank == 0):
      #print(strain)
      #print(total_force)
      #print("")
      graph_export.write(str(strain)+"  "+str(total_stress)+'\n')
         
    step += 1
#-------------------------------------------------------------------------------------------------------------------------------------


#-----------------------FINALITZATION OPERATIONS-------------------------------------------------------------------------------------- 

#

gid_io.FinalizeResults()
   
#os.chdir(graphs_path)
 
if(mpi.rank == 0): 
  graph_export.close() 

#os.chdir(list_path)

#
#

if(mpi.rank==0): 
  print 'Calculation ends at instant: ' + str(timer.time())
  elapsed_pr_time = timer.clock() - initial_pr_time
  elapsed_real_time = timer.time() - initial_real_time
  print 'Calculation ends at processing time instant: ' + str(timer.clock())
  print 'Elapsed processing time: ' + str(elapsed_pr_time)
  print 'Elapsed real time: ' + str(elapsed_real_time)
  print (my_timer)    
  print "COMPLETED ANALYSIS" 
 