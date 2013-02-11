# -*- coding: utf-8 -*-
import time as timer
import os
import sys
import math

import matplotlib
from numpy import *
from pylab import *  

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
#
#

from DEM_explicit_solver_var import *
from DEM_procedures import *


#defining a model part for the solid part

my_timer=Timer();
solid_model_part = ModelPart("SolidPart");  

import sphere_strategy as SolverStrategy
SolverStrategy.AddVariables(solid_model_part)

## reading the solid part: binary or ascii, multifile or single

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

#creating a solver object

solver = SolverStrategy.ExplicitStrategy(solid_model_part, domain_size) #here, solver variables initialize as default


#CREATING PATHS:

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

os.chdir(data_and_results)



if ( (ContinuumOption =="ON") and (ContactMeshOption =="ON") ) :
  
  contact_model_part = solver.contact_model_part   

  
ProcGiDSolverTransfer(solid_model_part,solver)

solver.Initialize()

dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

if(ModelDataInfo =="ON"):
  os.chdir(data_and_results)
  ProcModelData(solid_model_part,solver)       # calculates the mean number of neighbours the mean radius, etc..
  os.chdir(main_path)

if(ConcreteTestOption =="ON"):
  ProcListDefinition(solid_model_part,solver)  # defines the lists where we measure forces
  ProcSkinAndPressure(solid_model_part,solver)       # defines the skin and applies the pressure
  

if (CriticalTimeOption =="ON"):
  solver.Initial_Critical_Time() 

  if (dt!=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)):
    print("WARNING: Delta time has been modifyed to the critical one")
    dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

#initializations
time = 0.0
step = 0
time_old_print = 0.0

initial_pr_time = timer.clock()
initial_real_time = timer.time()

print('\n')
print ('Calculation starts at instant: ' + str(initial_pr_time)+'\n')

total_steps_expected = int(final_time/dt)
print ('Total number of TIME STEPs expected in the calculation are: ' + str(total_steps_expected) + ' if time step is kept ' +'\n' )



results = open('results.txt','w') #file to export some results
summary_results = open('summary_results.txt','w')

forcelist = []
forcelist2 = []
timelist = []
displacementlist = []

os.chdir(list_path)

multifile = open(problem_name+'_all'+'.post.lst','w')
multifile_5 = open(problem_name+'_5'+'.post.lst','w')
multifile_10 = open(problem_name+'_10'+'.post.lst','w')
multifile_50 = open(problem_name+'_50'+'.post.lst','w')


multifile.write('Multiple\n')
multifile_5.write('Multiple\n')
multifile_10.write('Multiple\n')
multifile_50.write('Multiple\n')

index_5 = 1
index_10 = 1
index_50 = 1

prev_time = 0.0
control = 0.0
cond = 0

os.chdir(main_path)

graph_export = open("strain_stress_data.csv",'w')
sigma_writting = open("mean_sigma.csv",'w')
sigma_writting2 = open("mean_sigma.csv",'w')

vs_radi = open("vs_radi.csv",'w')
vs_var_rad = open("vs_var_rad.csv",'w')



#Adding stress and strain lists
strainlist=[]
strainlist.append(0.0)
stresslist=[]
stresslist.append(0.0)

strain=0.0	

contact_model_part = solver.contact_model_part   

os.chdir(post_path)

if(Multifile == "single_file"):

  gid_io.InitializeMesh(0.0)

  gid_io.WriteMesh(contact_model_part.GetMesh());
  gid_io.FinalizeMesh()
  gid_io.InitializeResults(0.0, contact_model_part.GetMesh()); 

  gid_io.InitializeMesh(0.0)
  gid_io.WriteSphereMesh(solid_model_part.GetMesh())
  gid_io.FinalizeMesh()
  gid_io.InitializeResults(0.0, solid_model_part.GetMesh()); 

os.chdir(main_path)

# for the graph plotting    
velocity_node_y = 0.0
    
for node in sup_layer_fm:
    velocity_node_y = node.GetSolutionStepValue(VELOCITY_Y,0) #Applied velocity during the uniaxial compression test
    break
    
#done=False  #flag for the end of the confinement  
 
###################################################################
#                                                                 #
#--------------------------MAIN LOOP------------------------------#
#                                                                 #
###################################################################

while(time < final_time):
 
    dt = solid_model_part.ProcessInfo.GetValue(DELTA_TIME) #possible modifications of DELTA_TIME
    time = time + dt
    solid_model_part.CloneTimeStep(time)

    solid_model_part.ProcessInfo[TIME_STEPS] = step
        
    ####imprimint les forces en un arxiu.
    
    total_force=0
    force_node= 0
    
    os.chdir(data_and_results)
    
    for node in sup_layer_fm:
	
		force_node = node.GetSolutionStepValue(RHS,0)
		force_node_x = node.GetSolutionStepValue(RHS,0)[0]
		force_node_y = node.GetSolutionStepValue(RHS,0)[1]
		force_node_z = node.GetSolutionStepValue(RHS,0)[2]
	
	
		results.write(str(node.Id)+"  "+str(step)+"  "+str(force_node_y)+'\n')
		total_force += force_node_y


    #For a uniaxial compression test with a cylinder of 15 cm diameter and 30 cm height
    total_stress = total_force/(math.pi*75*75) #Stress in MPa
    stresslist.append(total_stress)

    #For a test tube of height 30 cm
    if(ContinuumOption =="ON"):
      strain += -2*velocity_node_y*dt/0.3
      strainlist.append(strain)
	   
    #writing lists to be printed
    forcelist.append(total_force)
    timelist.append(time)
    
    total_force=0
    force_node= 0
    
    for node in inf_layer_fm:
	
		force_node = node.GetSolutionStepValue(RHS,0)
		force_node_x = node.GetSolutionStepValue(RHS,0)[0]
		force_node_y = node.GetSolutionStepValue(RHS,0)[1]
		force_node_z = node.GetSolutionStepValue(RHS,0)[2]

		total_force += force_node_y
       
    #writing lists to be printed
    forcelist2.append(total_force)
        
    summary_results.write(str(step)+"  "+str(total_force)+'\n')

    os.chdir(main_path)
    
    #Dissable the confinement
    
    #if(ConcreteTestOption==True):
    
	  #if( (time > final_time*0.1) and (done==False)):
		#done=True;
		#for element in skin_list:
		  #element.SetValue(APPLIED_FORCE,(0,0,0))
		
		#print("Confinement finished at time "+str(time))
		
    solver.Solve()

   
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
	   	      

    os.chdir(list_path)
    
    multifile.write(problem_name+'_'+str(time)+'.post.bin\n')
    
    os.chdir(main_path)

  ##############     GiD IO        ################################################################################
    
    time_to_print = time - time_old_print
    #print str(time)
       
    if(time_to_print >= output_dt):
    
	os.chdir(graphs_path)

	#Drawing graph stress_strain:

	if( (ConcreteTestOption =="ON") and (RealTimeGraph =="ON") ):
	  clf()
	  plot(strainlist,stresslist,'b-')
	  grid(True)
	  title('Stress - Strain')
	  xlabel('Strain')
	  ylabel('Stress (MPa)')
	  savefig('Stress_strain') 

	os.chdir(main_path)	   	   

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

	  gid_io.InitializeMesh(time)
	  gid_io.WriteMesh(contact_model_part.GetMesh());
	  gid_io.FinalizeMesh()
	  gid_io.InitializeResults(time, contact_model_part.GetMesh()); 

	  
	##########PRINTING VARIABLES############
	
	ProcPrintingVariables(gid_io,solid_model_part,contact_model_part,time)  

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
	 
	index_5 += 1
	index_10 += 1
	index_50 += 1
	
	if(Multifile == "multiple_files"):
	  gid_io.FinalizeResults()
	
	os.chdir(main_path)     
              
	time_old_print = time
   
    #End of print loop

    os.chdir(main_path)
    
    graph_export.write(str(strain)+"  "+str(total_stress)+'\n')
      
    
    ####DEBUG ONLY # MIQUEL
    
    os.chdir(data_and_results)
    
    if(1<2):
    
      mean_sigma = solid_model_part.ProcessInfo[DOUBLE_DUMMY_2]
     
      sigma_writting2.write(str(time)+" "+str(mean_sigma)+'\n')
      
      if(total_stress > 1e-6):
          ratio_contact_total = mean_sigma/(1e6*total_stress);
          sigma_writting.write(str(time)+" "+str(ratio_contact_total)+'\n')
      
    #############
   
    step += 1

if(Multifile == "single_file"):
  gid_io.FinalizeResults()
  
  

os.chdir(data_and_results)

counter = 0
counter_all = 0

h   = 0.3
d   = 0.15
eps = 2


for element in solid_model_part.Elements:
  
  counter_all = counter_all +1

  node = element.GetNode(0)
  r = node.GetSolutionStepValue(RADIUS,0)
  x = node.X
  y = node.Y
  z = node.Z
  
  if ( ( (x*x+z*z) < ((d/2-eps*r)*(d/2-eps*r)) ) and ( (y>eps*r ) and (y<(h-eps*r)) ) ): 
    
    counter = counter + 1
    num_of_neigh = node.GetSolutionStepValue(NUM_OF_NEIGH,0)
    r = node.GetSolutionStepValue(RADIUS,0)
    volume_real = 3.141592*r*r*r*3*0.25
    volume_equ = node.GetSolutionStepValue(EQ_VOLUME_DEM,0)
    vs_radi.write(str(r)+"  "+str(volume_equ/volume_real)+'\n')
    vs_var_rad.write(str(num_of_neigh)+"  "+str(volume_equ/volume_real)+'\n')
  
  
#print("nodes interior")
#print (counter)

#print("nodes totals")
#print (counter_all)

#print("nodes skin")
#print (counter_all-counter)
  

total_volume = 0

for element in solid_model_part.Elements:
  
  node = element.GetNode(0)
  volume_equ = node.GetSolutionStepValue(REPRESENTATIVE_VOLUME,0)
  
  total_volume = total_volume + volume_equ
  
  

real_volume = 3.141592*d*d*0.25*h
  
print (total_volume)
print (real_volume)



print ( ( (total_volume/real_volume) - 1 ) * 100 )


vs_radi.close() 
vs_var_rad.close() 



graph_export.close() 
sigma_writting.close()
sigma_writting2.close()

results.close()
summary_results.close()

os.chdir(list_path)

multifile.close()
multifile_5.close()
multifile_10.close()
multifile_50.close()

###PLOTS

os.chdir(graphs_path)

if (1<2):
  clf()
  plot(timelist,forcelist,'b-')
  grid(True)
  title('Vertical force vs time')
  xlabel('time (s)')
  ylabel('Force (N)')
  #xlim(0.0,70000)
  #ylim(-5.0,103870403.214)
  #legend(('force'))
  savefig('Grafic_1')

if (1<2):
  clf()
  plot(timelist,forcelist2,'b-')
  grid(True)
  title('Vertical force vs time')
  xlabel('time (s)')
  ylabel('Force (N)')
  #xlim(0.0,70000)
  #ylim(-5.0,103870403.214)
  #legend(('force'))
  savefig('Grafic_2')
  
if (3<2):
  clf()
  plot(temps,Ep1,'b-', temps,Ec1,'r-',temps,Ee1,'g-',temps,Et1,'y-')
  grid(True)
  title('Energia potencial, cinetica, elastica i total en funcio del temps bola 1')
  xlabel('temps')
  ylabel('Energia')
  xlim(0.0,6.0)
  ylim(-5.0,250000.0)
  legend(('Ep','Ec','Ee','Et'))
  savefig('Grafic_energies_1_fe')

os.chdir(main_path)

print 'Calculation ends at instant: ' + str(timer.time())
elapsed_pr_time = timer.clock() - initial_pr_time
elapsed_real_time = timer.time() - initial_real_time
print 'Calculation ends at processing time instant: ' + str(timer.clock())
print 'Elapsed processing time: ' + str(elapsed_pr_time)
print 'Elapsed real time: ' + str(elapsed_real_time)
print (my_timer)    
print "COMPLETED ANALYSIS" 
