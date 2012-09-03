# -*- coding: utf-8 -*-
import DEM_explicit_solver_var
import time as timer

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
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions

gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_solid = ModelPartIO(input_file_name)
model_part_io_solid.ReadModelPart(solid_model_part)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
solid_model_part.SetBufferSize(2)

##adding dofs
SolverStrategy.AddDofs(solid_model_part)




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





#creating a solver object

dimension=DEM_explicit_solver_var.domain_size;

solver = SolverStrategy.ExplicitStrategy(solid_model_part, dimension); #here, solver variables initialize as default

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

solution_type = DEM_explicit_solver_var.SolutionType

if(solution_type == "Absolutal"):
    type_id = 2
else:
    type_id = 1

damp_ratio_type = DEM_explicit_solver_var.DampRatioType

if(damp_ratio_type == "ViscDamp"):
    damp_id = 2
elif(damp_ratio_type == "LocalDamp"):
    damp_id = 1
elif(damp_ratio_type == "BothDamp"):
    damp_id = 3
else:
    damp_id = 4
    
solver.damp_id=damp_id

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
rotation_spring_option	= DEM_explicit_solver_var.RotationalSpringOption

bounding_box_option 	= DEM_explicit_solver_var.BoundingBoxOption

global_variables_option = DEM_explicit_solver_var.GlobalVariablesOption
 

if(delta_option=="OFF"):
  search_radius_extension=0.0;

solver.time_scheme=time_scheme
solver.type_id=type_id

if(continuum_option =="ON"):
  solver.continuum_simulating_OPTION=True

if(delta_option =="ON"):
  solver.delta_OPTION=True
  
solver.search_radius_extension=search_radius_extension

if(rotation_option =="ON"):
  solver.rotation_OPTION=1  #xapuza
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

final_time = DEM_explicit_solver_var.max_time
output_dt  = DEM_explicit_solver_var.output_dt
dt = DEM_explicit_solver_var.max_time_step

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

dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

if (compute_critical_time =="ON"):
  solver.Critical_Time() 

  if (dt!=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)):
    print("WARNING: Delta time has been modifyed to the critical one")
    dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

#initializations
time = 0.0
step = 0
time_old_print = 0.0

current_pr_time = timer.clock()
current_real_time = timer.time()

print ('Calculation starts at instant: ' + str(current_pr_time)+'\n')

print ('Last TIME STEP is expected to be: ' + str(int(final_time/dt)) +'\n' )

results = open('results.txt','w') #file to export some results
summary_results = open('summary_results.txt','w')

forcelist = []
timelist = []
displacementlist = []

while(time < final_time):

    time = time + dt
    solid_model_part.CloneTimeStep(time)

    solid_model_part.ProcessInfo[TIME_STEPS] = step
    
    if (1<0):  ##this part is for a special test... will be erased.
      
	if (step > 4800):
	    print("python")
	    print(step)
	    print("len(Nodes)=",len(solid_model_part.Nodes))
	    a=None
	    count=0
	    for node in solid_model_part.Nodes:
	      if node.Id == 1890:
		a=node
		print "posicio:",count
	      count+=1
		
		#solid_model_part.Nodes[2896].GetSolutionStepValue(PARTICLE_FAILURE_ID,0)
	    print(a.GetSolutionStepValue(PARTICLE_FAILURE_ID,0)) 
	    print("python")
      
      ####imprimint les forces en un arxiu.
    
    total_force=0
    force_node= 0
    
    for node in force_measurement:
	
	force_node = node.GetSolutionStepValue(RHS,0)
	force_node_x = node.GetSolutionStepValue(RHS,0)[0]
	force_node_y = node.GetSolutionStepValue(RHS,0)[1]
	force_node_z = node.GetSolutionStepValue(RHS,0)[2]
	
	

      
	results.write(str(node.Id)+"  "+str(step)+"  "+str(force_node_y)+'\n')
	total_force += force_node_y
    
    
    #writing lists to be printed
    forcelist.append(total_force)
    timelist.append(time)
    
    
    summary_results.write(str(step)+"  "+str(total_force)+'\n')

       
    solver.Solve()

    #dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)
       
##############     GiD IO        ################################################################################
    time_to_print = time - time_old_print
    #print str(time)
    
    if(time_to_print >= DEM_explicit_solver_var.output_dt):
    
	print "TIME STEP = ", step
	gid_io.InitializeMesh(time);
        gid_io.WriteSphereMesh(solid_model_part.GetMesh());
        gid_io.FinalizeMesh();
	gid_io.InitializeResults(time, solid_model_part.GetMesh());   
        gid_io.WriteNodalResults(VELOCITY, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DISPLACEMENT, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(RHS, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(RADIUS, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PARTICLE_COHESION, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PARTICLE_TENSION, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PARTICLE_FAILURE_ID, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(GROUP_ID, solid_model_part.Nodes, time, 0)
        #gid_io.PrintOnGaussPoints(PARTICLE_FAILURE_ID, solid_model_part, time)   #there are no gauss points defined for spheres
        if (rotation_option == 1): ##xapuza
            gid_io.WriteNodalResults(ANGULAR_VELOCITY, solid_model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(MOMENT, solid_model_part.Nodes, time, 0)
        #gid_io.Flush()      
        gid_io.FinalizeResults()    
	time_old_print = time
    
   
    step += 1
    
results.close()
summary_results.close()


###PLOTS

import matplotlib
from numpy import *
from pylab import *  

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











print 'Calculation ends at instant: ' + str(timer.time())
elapsed_pr_time = timer.clock() - current_pr_time
elapsed_real_time = timer.time() - current_real_time
print 'Calculation ends at processing time instant: ' + str(timer.clock())
print 'Elapsed processing time: ' + str(elapsed_pr_time)
print 'Elapsed real time: ' + str(elapsed_real_time)
print (my_timer)    
print "COMPLETED ANALYSIS" 
