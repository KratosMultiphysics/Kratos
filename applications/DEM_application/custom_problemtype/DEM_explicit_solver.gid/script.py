# -*- coding: utf-8 -*-
import DEM_explicit_solver_var
import time as timer
import os
import sys
import math
import matplotlib
from numpy import *
from pylab import *  

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *


#defining a model part for the solid part
my_timer=Timer();
solid_model_part = ModelPart("SolidPart");  
#############################################

#introducing input file name
input_file_name = DEM_explicit_solver_var.problem_name
import sphere_strategy as SolverStrategy
SolverStrategy.AddVariables(solid_model_part)

#reading the solid part
gid_mode = GiDPostMode.GiD_PostBinary
#multifile = MultiFileFlag.MultipleFiles
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions

gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_solid = ModelPartIO(input_file_name)
model_part_io_solid.ReadModelPart(solid_model_part)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
solid_model_part.SetBufferSize(2)

##adding dofs
SolverStrategy.AddDofs(solid_model_part)

#creating a solver object

dimension=DEM_explicit_solver_var.domain_size

solver = SolverStrategy.ExplicitStrategy(solid_model_part, dimension) #here, solver variables initialize as default

##Obtaning options and values

integration_scheme = DEM_explicit_solver_var.Integration_Scheme
if (integration_scheme == 'forward_euler'):
    time_scheme = FowardEulerScheme()
elif (integration_scheme == 'mid_point_rule'):
    time_scheme = MidPointScheme()
elif (integration_scheme == 'const_average_acc'):
    time_scheme = ConstAverageAccelerationScheme()
else:
    print('scheme not defined')

normal_force_calculation = DEM_explicit_solver_var.NormalForceCalculation

if (normal_force_calculation == "Linear"):
    force_calculation_type_id = 0
elif (normal_force_calculation == "Hertz"):
    force_calculation_type_id = 1

damp_ID = DEM_explicit_solver_var.DampId

if(damp_ID == "ViscDamp"):
    damp_id = 1
else:
    damp_id = 0
    
solver.damp_id=damp_id

rota_damp_ID = DEM_explicit_solver_var.RotaDampId

if(rota_damp_ID == "LocalDamp"):
    rota_damp_id = 1
else:
    rota_damp_id = 0
    
solver.rota_damp_id=rota_damp_id


gravity = Vector(3)
gravity[0] = DEM_explicit_solver_var.gravity_x
gravity[1] = DEM_explicit_solver_var.gravity_y
gravity[2] = DEM_explicit_solver_var.gravity_z

solver.gravity=gravity

#options for the solver

compute_critical_time	= DEM_explicit_solver_var.CriticalTimeOption

continuum_option 	= DEM_explicit_solver_var.ContinuumOption
delta_option 		= DEM_explicit_solver_var.DeltaOption
search_radius_extension	= DEM_explicit_solver_var.search_radius_extension

rotation_option 	= DEM_explicit_solver_var.RotationOption
trihedron_option	= DEM_explicit_solver_var.TrihedronOption

rotation_spring_option	= DEM_explicit_solver_var.RotationalSpringOption

bounding_box_option 	= DEM_explicit_solver_var.BoundingBoxOption

global_variables_option = DEM_explicit_solver_var.GlobalVariablesOption

virtual_mass_option 	= DEM_explicit_solver_var.VirtualMassOption
virtual_mass_coeff 	= DEM_explicit_solver_var.VirtualMassCoefficient

contact_mesh_option = DEM_explicit_solver_var.ContactMeshOption
failure_criterion_option = DEM_explicit_solver_var.FailureCriterionOption

if(virtual_mass_option == "ON"):
    solver.virtual_mass_OPTION=1 #xapuza
    
solver.nodal_mass_coeff=virtual_mass_coeff    

if(delta_option=="OFF"):
  search_radius_extension=0.0;

solver.time_scheme=time_scheme
solver.force_calculation_type_id=force_calculation_type_id

if (compute_critical_time =="ON"):
  solver.critical_time_OPTION=1; #xapuza

if(delta_option =="ON"):
  solver.delta_OPTION=True

if(continuum_option =="ON"):
  solver.continuum_simulating_OPTION=True
  
  if(contact_mesh_option =="ON"):
    solver.contact_mesh_OPTION=1  #xapuza
  
  if(failure_criterion_option =="Mohr-Coulomb"):
    solver.failure_criterion_OPTION=1 
  elif(failure_criterion_option =="Uncoupled"):
    solver.failure_criterion_OPTION=2
    
  solver.tau_zero 		= DEM_explicit_solver_var.TauZero
  solver.sigma_max 		= DEM_explicit_solver_var.SigmaMax
  solver.sigma_min 		= DEM_explicit_solver_var.SigmaMin
  solver.internal_fricc 	= DEM_explicit_solver_var.InternalFricc
  
solver.search_radius_extension=search_radius_extension

if(rotation_option =="ON"):
  solver.rotation_OPTION=1  #xapuza
if(trihedron_option =="ON"):
  solver.trihedron_OPTION=1  #xapuza 
if(rotation_spring_option =="ON"):
  solver.rotation_spring_OPTION=1  #xapuza
       

       
solver.safety_factor = DEM_explicit_solver_var.dt_safety_factor #for critical time step calculation 

# global variable settings

if(global_variables_option =="ON"):
  solver.global_variables_OPTION = 1  #xapuza

solver.global_kn 	= DEM_explicit_solver_var.global_kn
solver.global_kt  	= DEM_explicit_solver_var.global_kt
solver.global_kr  	= DEM_explicit_solver_var.global_kr
solver.global_rn  	= DEM_explicit_solver_var.global_rn
solver.global_rt  	= DEM_explicit_solver_var.global_rt
solver.global_rr  	= DEM_explicit_solver_var.global_rr
solver.global_fri_ang 	= DEM_explicit_solver_var.global_fri_ang

# time settings

final_time 	= DEM_explicit_solver_var.max_time
output_dt  	= DEM_explicit_solver_var.output_dt
control_time	= DEM_explicit_solver_var.control_time
dt 		= DEM_explicit_solver_var.max_time_step

solver.delta_time=dt

# bounding box

n_step_destroy_distant = DEM_explicit_solver_var.search_step      # how many steps between elimination of distant particles?
n_step_search = DEM_explicit_solver_var.search_step
solver.n_step_search = n_step_search

bounding_box_enlargement_factor = DEM_explicit_solver_var.bounding_box_enlargement_factor    # by what factor do we want to enlarge the strict bounding box

extra_radius = 0.0
max_radius = 0.0
min_radius = 0.0
first_it = True

#calculation of search radius
for node in solid_model_part.Nodes:
      
  rad = node.GetSolutionStepValue(RADIUS)
  if rad > max_radius:  
      max_radius = rad
  if first_it == True:
      min_radius = rad
      first_it = False
  if rad < min_radius:  
      min_radius = rad
      
if(bounding_box_option =="ON"):
  solver.bounding_box_OPTION=1  #xapuza

extra_radius = 2.5 * max_radius
prox_tol = 0.000001 * min_radius  #currently not in use.
bounding_box_enlargement_factor = max(1.0 + extra_radius, bounding_box_enlargement_factor)

solver.enlargement_factor = bounding_box_enlargement_factor

#Initialize the problem.

solver.Initialize()

if(1<2):

  ## ADVANCED USER INTERACTION

  # Previous Calculations.

  Model_Data = open('Model_Data.txt','w')
  
  #mean radius, and standard deviation:

  i = 0
  sum_radi = 0
  sum_squared = 0
  for node in solid_model_part.Nodes:
  
    sum_radi += node.GetSolutionStepValue(RADIUS)
    sum_squared += node.GetSolutionStepValue(RADIUS)**2
    i+=1

  mean=sum_radi/i
  var =sum_squared/i-mean**2
  if(abs(var)<1e-05):
    var=0
  std_dev=var**0.5
  
  Model_Data.write("Radius Mean: "+str(mean)+'\n')
  Model_Data.write("Std Deviation: "+str(std_dev)+'\n')
  Model_Data.write('\n')
  
  Total_Particles 	= len(solid_model_part.Nodes)
  Total_Contacts  	= solver.model_part.ProcessInfo.GetValue(TOTAL_CONTACTS)/2
  Coordination_Number	= 1.0*(Total_Contacts*2)/Total_Particles
  
  Model_Data.write("Total Number of Particles: "+str(Total_Particles)+'\n')
  Model_Data.write("Total Number of Contacts: "+str(Total_Contacts)+'\n')
  Model_Data.write("Coordination Number NC: "+str(Coordination_Number)+'\n')
  Model_Data.write('\n')
  
  Model_Data.write("Volume Elements: "+str(DEM_explicit_solver_var.mass_elements)+'\n')
  
  Model_Data.close()

  # Defining lists (FOR COMPRESSION TESTS)

  sup_layer = list()
  inf_layer = list()
  fix_particles = list()
  force_measurement = list()
  special_selection = list()
  others = list()
  
  for node in solid_model_part.Nodes:
    if (node.GetSolutionStepValue(GROUP_ID)==1):
      sup_layer.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==2):
      inf_layer.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==3):
      fix_particles.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==4):
      force_measurement.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==5):
      special_selection.append(node)
    else:
      others.append(node)


dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

if (compute_critical_time =="ON"):
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


#paths:

main_path 	 = os.getcwd()
post_path 	 = str(main_path)+'/'+str(input_file_name)+'_Post_Files'
list_path 	 = str(main_path)+'/'+str(input_file_name)+'_Post_Lists'
neigh_list_path  = str(main_path)+'/'+str(input_file_name)+'_Neigh_Lists'
data_and_results = str(main_path)+'/'+str(input_file_name)+'_Results_and_Data'
graphs_path	 = str(main_path)+'/'+str(input_file_name)+'_Graphs'	
MPI_results    = str(main_path)+'/'+str(input_file_name)+'_MPI_results'	

for directory in [post_path, list_path, neigh_list_path, data_and_results, graphs_path, MPI_results]:

  if not os.path.isdir(directory):
    
      os.makedirs(str(directory))

os.chdir(data_and_results)

results = open('results.txt','w') #file to export some results
summary_results = open('summary_results.txt','w')

forcelist = []
forcelist2 = []
timelist = []
displacementlist = []

os.chdir(list_path)

multifile = open(input_file_name+'_all'+'.post.lst','w')
multifile_5 = open(input_file_name+'_5'+'.post.lst','w')
multifile_10 = open(input_file_name+'_10'+'.post.lst','w')
multifile_50 = open(input_file_name+'_50'+'.post.lst','w')


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

#Adding stress and strain lists
strainlist=[]
strainlist.append(0.0)
stresslist=[]
stresslist.append(0.0)

strain=0.0	

contact_model_part = solver.contact_model_part   

os.chdir(post_path)

gid_io.InitializeMesh(0.0)
gid_io.WriteSphereMesh(solid_model_part.GetMesh())
gid_io.FinalizeMesh()
gid_io.InitializeResults(0.0, solid_model_part.GetMesh()); 

gid_io.InitializeMesh(0.0)
gid_io.WriteMesh(contact_model_part.GetMesh());
gid_io.FinalizeMesh()
gid_io.InitializeResults(0.0, contact_model_part.GetMesh()); 

os.chdir(main_path)

#Defining list of skin particles (For a test tube of height 30 cm and diameter 15 cm)
    
for element in solid_model_part.Elements:
	
    element.SetValue(SKIN_SPHERE,0)
    node = element.GetNode(0)
    x = node.X
    y = node.Y
    z = node.Z
    r = node.GetSolutionStepValue(RADIUS,0)
    h=0.3
    d=0.15
    eps=2
	
    if ( (x*x+z*z)>=((d/2-eps*r)*(d/2-eps*r)) or y<=eps*r or y>=(h-eps*r) ): #For a cylinder test with the center of the base at (0,0,0) and with the axis y normal to the base
	element.SetValue(SKIN_SPHERE,1)
    
    velocity_node_y = 0.0
    
for node in force_measurement:
    velocity_node_y = node.GetSolutionStepValue(VELOCITY_Y,0) #Applied velocity during the uniaxial compression test

while(time < final_time):
  
  
    dt = solid_model_part.ProcessInfo.GetValue(DELTA_TIME) #possible modifications of DELTA_TIME
    time = time + dt
    solid_model_part.CloneTimeStep(time)

    solid_model_part.ProcessInfo[TIME_STEPS] = step
        
    ####imprimint les forces en un arxiu.
    
    total_force=0
    force_node= 0
    
    os.chdir(data_and_results)
    
    for node in force_measurement:
	
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
    if(continuum_option =="ON"):
      strain += -velocity_node_y*dt/0.3
      strainlist.append(strain)

	   
    #writing lists to be printed
    forcelist.append(total_force)
    timelist.append(time)
    
    total_force=0
    force_node= 0
    
    for node in inf_layer:
	
	force_node = node.GetSolutionStepValue(RHS,0)
	force_node_x = node.GetSolutionStepValue(RHS,0)[0]
	force_node_y = node.GetSolutionStepValue(RHS,0)[1]
	force_node_z = node.GetSolutionStepValue(RHS,0)[2]

	total_force += force_node_y
    
    
    #writing lists to be printed
    forcelist2.append(total_force)
    
    
    summary_results.write(str(step)+"  "+str(total_force)+'\n')

    os.chdir(main_path)
    
    solver.Solve()

    #dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)
    
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
    
    multifile.write(input_file_name+'_'+str(time)+'.post.bin\n')
    
    os.chdir(main_path)
    
##############     GiD IO        ################################################################################
    time_to_print = time - time_old_print
    #print str(time)
    
    if(time_to_print >= DEM_explicit_solver_var.output_dt):
    
	os.chdir(graphs_path)

	#Drawing graph stress_strain:

	clf()
	plot(strainlist,stresslist,'b-')
	grid(True)
	title('Stress - Strain')
	xlabel('Strain')
	ylabel('Stress (MPa)')
	savefig('Stress_strain') 

	os.chdir(main_path)	   	   

	if(2<3): #printing neighbours id's
	  
   
	  os.chdir(neigh_list_path)
	  neighbours_list = open('neigh_list_'+ str(time),'w')
      
	  for elem in solid_model_part.Elements:
	
	      ID=(elem.Id)
	      #Neigh_ID = elem.GetValue(NEIGHBOURS_IDS_DOUBLE)
	      Neigh_ID = elem.GetValue(NEIGHBOURS_IDS)
	      #print(len(Neigh_ID))
	      
	      for i in range(len(Neigh_ID)):
	
		neighbours_list.write(str(ID)+' '+str(Neigh_ID[i])+'\n')
	  
	  neighbours_list.close()
	  os.chdir(main_path)
	
	os.chdir(post_path)
	
	#print "TIME STEP = ", step
	#gid_io.InitializeMesh(time);
        #gid_io.WriteSphereMesh(solid_model_part.GetMesh());
        #gid_io.FinalizeMesh();
	#gid_io.InitializeResults(time, solid_model_part.GetMesh());
	
	
	if (DEM_explicit_solver_var.print_velocity=="1"):
	  gid_io.WriteNodalResults(VELOCITY, contact_model_part.Nodes, time, 0)	  
        if (DEM_explicit_solver_var.print_displacement=="1"):
	  gid_io.WriteNodalResults(DISPLACEMENT, contact_model_part.Nodes, time, 0)       
        if (DEM_explicit_solver_var.print_rhs=="1"):
	  gid_io.WriteNodalResults(RHS, contact_model_part.Nodes, time, 0)       
        if (DEM_explicit_solver_var.print_total_forces=="1"):
	  gid_io.WriteNodalResults(TOTAL_FORCES, contact_model_part.Nodes, time, 0)	  
        if (DEM_explicit_solver_var.print_damp_forces=="1"):
	  gid_io.WriteNodalResults(DAMP_FORCES, contact_model_part.Nodes, time, 0)        
        if (DEM_explicit_solver_var.print_radius=="1"):
	  gid_io.WriteNodalResults(RADIUS, contact_model_part.Nodes, time, 0)       
        if (DEM_explicit_solver_var.print_particle_cohesion=="1"):
	  gid_io.WriteNodalResults(PARTICLE_COHESION, contact_model_part.Nodes, time, 0)       
        if (DEM_explicit_solver_var.print_particle_tension=="1"):
	  gid_io.WriteNodalResults(PARTICLE_TENSION, contact_model_part.Nodes, time, 0)
        if (DEM_explicit_solver_var.print_group_id=="1"):
	  gid_io.WriteNodalResults(GROUP_ID, contact_model_part.Nodes, time, 0)
        if (DEM_explicit_solver_var.print_export_id=="1"):
	  gid_io.WriteNodalResults(EXPORT_ID, contact_model_part.Nodes, time, 0)
        if (DEM_explicit_solver_var.print_export_particle_failure_id=="1"):
	  gid_io.WriteNodalResults(EXPORT_PARTICLE_FAILURE_ID, contact_model_part.Nodes, time, 0)
	if (DEM_explicit_solver_var.print_export_skin_sphere=="1"):
	  gid_io.WriteNodalResults(EXPORT_SKIN_SPHERE, contact_model_part.Nodes, time, 0)
	
	if (contact_mesh_option == "ON"): ##xapuza
	  if (DEM_explicit_solver_var.print_local_contact_force_low=="1"):
	    gid_io.PrintOnGaussPoints(LOCAL_CONTACT_FORCE_LOW,contact_model_part,time)
	  if (DEM_explicit_solver_var.print_local_contact_force_high=="1"):
	    gid_io.PrintOnGaussPoints(LOCAL_CONTACT_FORCE_HIGH,contact_model_part,time)
	  if (DEM_explicit_solver_var.print_contact_failure=="1"): 
	    gid_io.PrintOnGaussPoints(CONTACT_FAILURE,contact_model_part,time)	 
	  if (DEM_explicit_solver_var.print_failure_criterion_state=="1"):
	    gid_io.PrintOnGaussPoints(FAILURE_CRITERION_STATE,contact_model_part,time)  	    
	  if (DEM_explicit_solver_var.print_contact_tau=="1"):
	    gid_io.PrintOnGaussPoints(CONTACT_TAU,contact_model_part,time)
	  if (DEM_explicit_solver_var.print_contact_sigma=="1"):
	    gid_io.PrintOnGaussPoints(CONTACT_SIGMA,contact_model_part,time)
	  
             
        if (rotation_option == "ON"): ##xapuza
            if (DEM_explicit_solver_var.print_angular_velocity=="1"):
	      gid_io.WriteNodalResults(ANGULAR_VELOCITY, contact_model_part.Nodes, time, 0)
            if (DEM_explicit_solver_var.print_particle_moment=="1"):
	      gid_io.WriteNodalResults(PARTICLE_MOMENT, contact_model_part.Nodes, time, 0)
           # if (DEM_explicit_solver_var.print_euler_angles):
	      #gid_io.WriteLocalAxesOnNodes(EULER_ANGLES, contact_model_part.Nodes, time, 0)
        
        #gid_io.FinalizeResults() 
        gid_io.Flush()
        sys.stdout.flush()
        
        os.chdir(data_and_results)
        
        if (index_5==5):
	  
	  multifile_5.write(input_file_name+'_'+str(time)+'.post.bin\n')
	  
	  index_5=0
	  
	if (index_10==10):
	  
	  multifile_10.write(input_file_name+'_'+str(time)+'.post.bin\n')
	  
	  index_10=0
	  
	if (index_50==50):
	  
	  multifile_50.write(input_file_name+'_'+str(time)+'.post.bin\n')
	  
	  index_50=0
	 
	index_5 += 1
	index_10 += 1
	index_50 += 1
	
	os.chdir(main_path)
       
              
           
	time_old_print = time
    

    os.chdir(main_path)
    
    graph_export.write(str(strain)+"  "+str(total_stress)+'\n')
    
    step += 1

  
    
gid_io.FinalizeResults()

os.chdir(data_and_results)

graph_export.close() 
results.close()
summary_results.close()

os.chdir(list_path)

multifile.close()
multifile_5.close()
multifile_10.close()
multifile_50.close()

os.chdir(graphs_path)

###PLOTS


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

clf()
plot(strainlist,stresslist,'b-')
grid(True)
title('Stress - Strain')
xlabel('Strain')
ylabel('Stress (MPa)')
savefig('Stress_strain')
  
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
