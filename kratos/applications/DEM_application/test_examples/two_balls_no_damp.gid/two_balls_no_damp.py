import time as timer
import os
import sys
import math

import matplotlib
from numpy import *
from pylab import *  

##import cProfile
##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_path = '../../../..'         #########################NUEVO
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
kratos_benchmarking_path = '../../../../benchmarking' #######NUEVO

import sys
sys.path.append(kratos_path)
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_benchmarking_path) ###################NUEVO

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

from DEM_explicit_solver_var import *
from DEM_procedures import *

import benchmarking                 #########################NUEVO

def FindNode(node_list,X,Y,Z):
    for node in node_list:
	if((node.X-X)**2 + (node.Y-Y)**2 + (node.Z-Z)**2 < .000001):
		return node
    
def BenchmarkCheck(time, node1):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_Y), "Node Displacement", 1.0)

#---------------------MODEL PART KRATOS AND GID.IO ------------------------------------------------------------------

#defining a model part for the solid part
my_timer=Timer();
solid_model_part = ModelPart("SolidPart");  

import sphere_strategy as SolverStrategy
SolverStrategy.AddVariables(solid_model_part)

## reading the solid part: binary or ascii, multifile or single --> only binary and single for mpi.

if(OutputFileType == "Binary"):
  gid_mode = GiDPostMode.GiD_PostBinary
else:
  gid_mode = GiDPostMode.GiD_PostAscii
  
if(Multifile == "multiple_files"):
  multifile = MultiFileFlag.MultipleFiles
else:
  multifile = MultiFileFlag.SingleFile

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions

gid_io = GidIO(problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_solid = ModelPartIO(problem_name)

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

main_path 	 = os.getcwd()
post_path 	 = str(main_path)+'/'+str(problem_name)+'_Post_Files'
list_path 	 = str(main_path)+'/'+str(problem_name)+'_Post_Lists'
neigh_list_path  = str(main_path)+'/'+str(problem_name)+'_Neigh_Lists'
data_and_results = str(main_path)+'/'+str(problem_name)+'_Results_and_Data'
graphs_path	 = str(main_path)+'/'+str(problem_name)+'_Graphs'	
MPI_results    = str(main_path)+'/'+str(problem_name)+'_MPI_results'	

for directory in [post_path, list_path, neigh_list_path, data_and_results, graphs_path, MPI_results]:

  if not os.path.isdir(directory):
    
      os.makedirs(str(directory))

os.chdir(list_path)

multifile = open(problem_name+'_all'+'.post.lst','w'); multifile_5 = open(problem_name+'_5'+'.post.lst','w');
multifile_10 = open(problem_name+'_10'+'.post.lst','w'); multifile_50 = open(problem_name+'_50'+'.post.lst','w')

multifile.write('Multiple\n'), multifile_5.write('Multiple\n'); multifile_10.write('Multiple\n'); multifile_50.write('Multiple\n')
index_5 = 1; index_10 = 1; index_50 = 1; prev_time = 0.0; control = 0.0; cond = 0

os.chdir(main_path)

graph_export = open("strain_stress_data.csv",'w');

#Adding stress and strain lists
strainlist=[]; strainlist.append(0.0)
stresslist=[]; stresslist.append(0.0)
strain=0.0; total_stress = 0.0; first_time_entry = 1

# for the graph plotting    
velocity_node_y = 0.0
    
for node in sup_layer_fm:
    velocity_node_y = node.GetSolutionStepValue(VELOCITY_Y,0) #Applied velocity during the uniaxial compression test
    print 'velocity for the graph' + str(velocity_node_y) + '\n'
    break

export_model_part = solid_model_part

if ( (ContinuumOption =="ON") and (ContactMeshOption =="ON") ) :
  
  contact_model_part = solver.contact_model_part   
  export_model_part = contact_model_part
 
#-------------------------------------------------------------------------------------------------------------------------

#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALITZATION--------------------------------------------------------

Pressure = 0.0  
Pressure = ProcGiDSolverTransfer(solid_model_part,solver)

if(ModelDataInfo =="ON"):
  os.chdir(data_and_results)
  ProcModelData(solid_model_part,solver)       # calculates the mean number of neighbours the mean radius, etc..
  os.chdir(main_path)

if(ConcreteTestOption =="ON"):
  ProcListDefinition(solid_model_part,solver)  # defines the lists where we measure forces
  (SKIN, LAT, BOT, TOP, XLAT, XTOP, XBOT, XTOPCORNER, XBOTCORNER) = ProcSkinAndPressure(solid_model_part,solver)       # defines the skin and applies the pressure

#mesurement
heigh = 0.3

if(ContinuumOption =="ON" and ConcreteTestOption =="ON"):
  
  Y_mean_bot = ProcMeasureBOT(BOT,solver)
  Y_mean_top = ProcMeasureTOP(TOP,solver)
  ini_heigh = Y_mean_top - Y_mean_bot
  
  print ('Initial Heigh of the Model: ' + str(ini_heigh)+'\n')
  heigh = ini_heigh

print'Initialitzating Problem....'
solver.Initialize()
print 'Initialitzation Complete' + '\n'

###############################################################
node1 = FindNode(solid_model_part.Nodes , 0.0, 1.0, 0.0)
print node1 #there is a memory problem with the string
###############################################################

dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

if (CriticalTimeOption =="ON"):
  solver.Initial_Critical_Time() 

  if (dt!=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)):
    print("WARNING: Delta time has been modifyed to the critical one")
    dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

time = 0.0; step = 0; time_old_print = 0.0

initial_pr_time = timer.clock()
initial_real_time = timer.time()

print ('SOLVE starts at instant: ' + str(initial_pr_time)+'\n')

total_steps_expected = int(final_time/dt)
print ('Total number of TIME STEPs expected in the calculation are: ' + str(total_steps_expected) + ' if time step is kept ' +'\n' )
    
#-------------------------------------------------------------------------------------------------------------------------------------

#-----------------------SINGLE FILE MESH AND RESULTS INITIALITZATION-------------------------------------------------------------------


os.chdir(post_path)

if(Multifile == "single_file"):

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

os.chdir(main_path)
while(time < final_time):
 
    dt = solid_model_part.ProcessInfo.GetValue(DELTA_TIME) #possible modifications of DELTA_TIME
    time = time + dt
    solid_model_part.CloneTimeStep(time)

    solid_model_part.ProcessInfo[TIME_STEPS] = step

    #########################_SOLVE_#########################################4
    
    os.chdir(main_path) 
    solver.Solve()
    
    #########################TIME CONTROL######################################4
   
    incremental_time = (timer.time()-initial_real_time)- prev_time

    if (incremental_time > control_time): 
      
      percentage = 100.0*(float(step)/total_steps_expected)
      print 'Real time calculation: ' + str(timer.time()-initial_real_time) 
      print 'Percentage Completed: ' +str(percentage) + ' %' 
      print "TIME STEP = " + str(step) + '\n'
        
      prev_time = (timer.time()-initial_real_time)
    
    if ( (timer.time()-initial_real_time > 60.0) and cond==0):
    
      cond=1
    
      estimation_time=60.0*(total_steps_expected/step) #seconds
    
      print('the total calculation estimated time is '+str(estimation_time)+'seconds.'+'\n')
      print('in minutes :'+str(estimation_time/60)+'min.'+'\n')
      print('in hours :'+str((estimation_time/60)/60)+'hrs.'+'\n')
      print('in days :'+str(((estimation_time/60)/60)/24)+'days.'+'\n') 
    
      if (((estimation_time/60)/60)/24 > 2.0):
        print('WARNING!!!:       VERY LASTING CALCULATION'+'\n')
    
    #########################CONCRETE_TEST_STUFF#########################################4
 
    if( (ConcreteTestOption =="ON") and (step==3) ):
      
      #Cross section Area Control
      print '\n' + '----------------------CONCRETE TEST CONTROLS----------------------' + '\n'
      print 'Total Horitzontal Numerical Cross Section on Force Measurement: ' + str(solid_model_part.ProcessInfo.GetValue(AREA_VERTICAL_TAPA))
      #print( solid_model_part.ProcessInfo.GetValue(AREA_VERTICAL_CENTRE) )
      
      total_volume = 0.0;  h   = 0.3;    d   = 0.15
      
      for element in solid_model_part.Elements:
  
        node = element.GetNode(0)
        volume_equ = node.GetSolutionStepValue(REPRESENTATIVE_VOLUME,0) 
        total_volume += volume_equ

      real_volume = 3.141592*d*d*0.25*h
      
      print 'Total Numerical Volume: ' + str(total_volume)
      print 'Total Numerical Volume: ' + str(real_volume)
      print 'Error: ' + str(100*abs(total_volume-real_volume)/real_volume) +'%'+'\n'
      print '------------------------------------------------------------------' + '\n'
    
    os.chdir(data_and_results)
    
    total_force=0
    force_node= 0
    
    #For a uniaxial compression test with a cylinder of 15 cm diameter and 30 cm height

    if( ContinuumOption =="ON" and ( time > 0.01*TimePercentageFixVelocities*final_time) and ConcreteTestOption =="ON" and ConcreteTestOption =="ON" ):
    
      if(first_time_entry):
        Y_mean_bot = ProcMeasureBOT(BOT,solver)
        Y_mean_top = ProcMeasureTOP(TOP,solver)
        
        ini_heigh2 = Y_mean_top - Y_mean_bot
        
        print 'Current Heigh after confinement: ' + str(ini_heigh2) + '\n'
        print 'Axial strain due to the confinement: ' + str( 100*(ini_heigh2-ini_heigh)/ini_heigh ) + ' %' +'\n'

        first_time_entry = 0
        heigh = ini_heigh2
        
      strain += -2*velocity_node_y*dt/heigh
      strainlist.append(strain)
      
      for node in sup_layer_fm:
      
        force_node = node.GetSolutionStepValue(RHS,0)
        force_node_x = node.GetSolutionStepValue(RHS,0)[0]
        force_node_y = node.GetSolutionStepValue(RHS,0)[1]
        force_node_z = node.GetSolutionStepValue(RHS,0)[2]
        
        total_force += force_node_y
      
      total_stress = total_force/(math.pi*75*75) + (1e-6)*Pressure #Stress in MPa
      stresslist.append(total_stress)

    os.chdir(list_path)
    
    multifile.write(problem_name+'_'+str(time)+'.post.bin\n')
    
    os.chdir(main_path)

  #########################___GiD IO____#########################################4

    time_to_print = time - time_old_print

    if(time_to_print >= output_dt):

      BenchmarkCheck(time, node1)
      os.chdir(data_and_results)
        
      if (index_5==5):
        
        multifile_5.write(problem_name+'_'+str(time)+'.post.bin\n')
        index_5=0
        
      if (index_10==10):
        
        multifile_10.write(problem_name+'_'+str(time)+'.post.bin\n')
        index_10=0
        
      if (index_50==50):
        multifile_50.write(problem_name+'_'+str(time)+'.post.bin\n')
        index_50=0
      
      index_5 += 1;      index_10 += 1;      index_50 += 1

      if(Multifile == "multiple_files"):
        gid_io.FinalizeResults()

      os.chdir(graphs_path)

      #Drawing graph stress_strain:

      if( (ConcreteTestOption =="ON") and (RealTimeGraph =="ON") ):
        
        clf()
        plot(strainlist,stresslist,'b-')
        grid(True)
        title('Stress - Strain')
        xlabel('Total Vertical Strain')
        ylabel('Total Stress (MPa)')
        savefig('Stress_strain') 

      if(PrintNeighbourLists == "ON"): #printing neighbours id's
  
        os.chdir(neigh_list_path)
        neighbours_list = open('neigh_list_'+ str(time),'w')
  
        for elem in solid_model_part.Elements:
          ID=(elem.Id)
          Neigh_ID = elem.GetValue(NEIGHBOURS_IDS)
      
          for i in range(len(Neigh_ID)):

            neighbours_list.write(str(ID)+' '+str(Neigh_ID[i])+'\n')
  
        neighbours_list.close()

      os.chdir(post_path)

      if(Multifile == "multiple_files"):

        gid_io.InitializeMesh(time)
        gid_io.WriteSphereMesh(solid_model_part.GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(time, solid_model_part.GetMesh()); 
      
        if (ContactMeshOption =="ON"): 
          gid_io.InitializeMesh(time)
          gid_io.WriteMesh(contact_model_part.GetMesh());
          gid_io.FinalizeMesh()
          gid_io.InitializeResults(time, contact_model_part.GetMesh()); 
 
      ProcPrintingVariables(gid_io,export_model_part,time)  

      os.chdir(main_path)     
              
      time_old_print = time

    graph_export.write(str(strain)+"  "+str(total_stress)+'\n')
         
    step += 1
#-------------------------------------------------------------------------------------------------------------------------------------


#-----------------------FINALITZATION OPERATIONS-------------------------------------------------------------------------------------- 

if(Multifile == "single_file"):

  gid_io.FinalizeResults()
   
os.chdir(graphs_path)
 
graph_export.close() 

os.chdir(list_path)

multifile.close(); multifile_5.close(); multifile_10.close(); multifile_50.close()
os.chdir(main_path)

print 'Calculation ends at instant: ' + str(timer.time())
elapsed_pr_time = timer.clock() - initial_pr_time
elapsed_real_time = timer.time() - initial_real_time
print 'Calculation ends at processing time instant: ' + str(timer.clock())
print 'Elapsed processing time: ' + str(elapsed_pr_time)
print 'Elapsed real time: ' + str(elapsed_real_time)
print (my_timer)    
print "COMPLETED ANALYSIS" 
