import fluid_ulf_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = fluid_ulf_var.domain_size

##################################################################
##################################################################
# ATTENTION: here the order is important

import time as timer
import os
import sys
import math
import matplotlib
from numpy import *
from pylab import * 

#including kratos path
sys.path.append(fluid_ulf_var.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.PFEMApplication import PfemUtils
from KratosMultiphysics.StructuralApplication import *
##############################################################################################                   
#   (_)(_)(_)(_)            (_)(_)(_)(_)(_)         (_) _     _ (_)   
#    (_)      (_)_          (_)                     (_)(_)   (_)(_)   
#    (_)        (_)         (_) _  _                (_) (_)_(_) (_)   
#    (_)        (_)         (_)(_)(_)               (_)   (_)   (_)   
#    (_)       _(_)         (_)                     (_)         (_)   
#    (_)_  _  (_)           (_) _  _  _  _          (_)         (_)   
#   (_)(_)(_)(_)            (_)(_)(_)(_)(_)         (_)         (_)   
################################################################################################                                                                                                                                                           

# -*- coding: utf-8 -*-
from KratosMultiphysics.DEMApplication import *
from DEM_explicit_solver_var import *
from DEM_procedures import *                                                                                                                                                          
################################################################################################
#                #
#               ###                
#              #####              
#             #######
#            #########
#               ###
#               ###     
#################################################################################################
#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");
structure_model_part = ModelPart("StructurePart");  
combined_model_part = ModelPart("CombinedPart");

SolverType=fluid_ulf_var.SolverType
if (SolverType=="Incompressible_Modified_FracStep" or SolverType=="FracStep"):
    fluid_only_model_part = ModelPart("FluidOnlyPart");

#############################################
#importing the solvers needed
if(SolverType == "Incompressible_Modified_FracStep"):
    import ulf_frac
    ulf_frac.AddVariables(fluid_model_part)   
elif(SolverType == "FracStep"):
    import ulf_frac
    ulf_frac.AddVariables(fluid_model_part)       
elif(SolverType == "Quasi_Inc_Constant_Pressure"):
    import ulf_fsi
    ulf_fsi.AddVariables(fluid_model_part)
elif(SolverType == "Quasi_Inc_Linear_Pressure"):
    import ulf_fsi_inc
    ulf_fsi_inc.AddVariables(fluid_model_part)
else:
    raise "solver type not supported: options are fractional_step - \
	modified_frac_steop - quasi_inc_constant_pres - \
	quasi_inc_lin_pres"

#introducing input file name
input_file_name = fluid_ulf_var.problem_name

##############################################################################################                   
#   (_)(_)(_)(_)            (_)(_)(_)(_)(_)         (_) _     _ (_)   
#    (_)      (_)_          (_)                     (_)(_)   (_)(_)   
#    (_)        (_)         (_) _  _                (_) (_)_(_) (_)   
#    (_)        (_)         (_)(_)(_)               (_)   (_)   (_)   
#    (_)       _(_)         (_)                     (_)         (_)   
#    (_)_  _  (_)           (_) _  _  _  _          (_)         (_)   
#   (_)(_)(_)(_)            (_)(_)(_)(_)(_)         (_)         (_)   
################################################################################################   
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
  DEM_multifile = MultiFileFlag.MultipleFiles
else:
  DEM_multifile = MultiFileFlag.SingleFile

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
DEM_gid_io = GidIO(problem_name, gid_mode, DEM_multifile, deformed_mesh_flag, write_conditions)
#model_part_io_solid = ModelPartIO(problem_name)
#model_part_io_solid.ReadModelPart(solid_model_part)
################################################################################################
#                #
#               ###                
#              #####              
#             #######
#            #########
#               ###
#               ###     
#################################################################################################
#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
#model_part_io_origin = ModelPartIO(input_file_name)
#model_part_io_origin.ReadModelPart(fluid_model_part)
#################################################################################################    _   __            
#   / | / /__ _      __
#  /  |/ / _ \ | /| / /
# / /|  /  __/ |/ |/ / 
#/_/ |_/\___/|__/|__/ 
#################################################################################################
#reading the fluid model part
model_part_io_origin = ModelPartIO(input_file_name)
model_part_io_origin.ReadModelPart(fluid_model_part)
print fluid_model_part
print "ULF model read correctly"

#reading the DEM model part
data_io = ModelPartIO(problem_name)
print problem_name
data_io.ReadModelPart(solid_model_part)
print solid_model_part
print "DEM model read correctly"

################################################################################################# 
#| |     / /__  ____   
#| | /| / / _ \/ __ \  
#| |/ |/ /  __/ / / /  
#|__/|__/\___/_/ /_/   
#################################################################################################                     
#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

compute_reactions=fluid_ulf_var.compute_reactions
##adding dofs
if(SolverType == "Incompressible_Modified_FracStep"):
    ulf_frac.AddDofs(fluid_model_part, compute_reactions)
elif(SolverType == "FracStep"):
    ulf_frac.AddDofs(fluid_model_part, compute_reactions)
elif(SolverType == "Quasi_Inc_Constant_Pressure"):
    ulf_fsi.AddDofs(fluid_model_part, compute_reactions)
elif(SolverType == "Quasi_Inc_Linear_Pressure"):
    ulf_fsi_inc.AddDofs(fluid_model_part, compute_reactions)

if(SolverType == "Quasi_Inc_Constant_Pressure" or SolverType == "Quasi_Inc_Linear_Pressure"):
      for node in fluid_model_part.Nodes:
	  node.Free(PRESSURE)

#setting the limits of the bounding box
box_corner1 = Vector(3); 
box_corner1[0]=fluid_ulf_var.bounding_box_corner1_x; box_corner1[1]=fluid_ulf_var.bounding_box_corner1_y; box_corner1[2]=fluid_ulf_var.bounding_box_corner1_z;
box_corner2 = Vector(3); 
box_corner2[0]=fluid_ulf_var.bounding_box_corner2_x; box_corner2[1]=fluid_ulf_var.bounding_box_corner2_y; box_corner2[2]=fluid_ulf_var.bounding_box_corner2_z;

#here we write the convergence data..,
outstring2 = "convergence_info.txt"
outputfile1 = open(outstring2, 'w')

add_nodes=fluid_ulf_var.adaptive_refinement
bulk_modulus=fluid_ulf_var.bulk_modulus
density=fluid_ulf_var.density
FSI=fluid_ulf_var.FSI
#creating the solvers
#fluid solver
##check to ensure that no node has zero density or pressure
is_fsi_interf=0.0

if(SolverType == "Incompressible_Modified_FracStep"):    
    solver = ulf_frac.ULF_FSISolver(outputfile1, fluid_only_model_part, fluid_model_part, structure_model_part, combined_model_part, FSI, compute_reactions, box_corner1, box_corner2, domain_size, add_nodes, bulk_modulus, density)
    solver.alpha_shape = fluid_ulf_var.alpha_shape;
    solver.echo_level = 2;
    
    for node in fluid_model_part.Nodes:
	node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)
	node.SetSolutionStepValue(DENSITY,0, density)
	node.SetSolutionStepValue(VISCOSITY,0, 0.0001)
	node.SetSolutionStepValue(BODY_FORCE_Y,0, -10.0001) 
    solver.Initialize()
    
if(SolverType == "FracStep"):    
    solver = ulf_frac.ULF_FSISolver(outputfile1, fluid_only_model_part, fluid_model_part, structure_model_part, combined_model_part, FSI, compute_reactions, box_corner1, box_corner2, domain_size, add_nodes, bulk_modulus, density)
    solver.alpha_shape = fluid_ulf_var.alpha_shape;
    solver.echo_level = 2;
    for node in fluid_model_part.Nodes:
	node.SetSolutionStepValue(BULK_MODULUS,0, 0.0)   
	node.SetSolutionStepValue(DENSITY,0, density)   
    solver.Initialize()
    
elif(SolverType == "Quasi_Inc_Constant_Pressure"):
    solver = ulf_fsi.ULF_FSISolver(fluid_model_part, structure_model_part, combined_model_part, compute_reactions, box_corner1, box_corner2, domain_size, add_nodes)
    solver.alpha_shape = fluid_ulf_var.alpha_shape;
    solver.echo_level = 2;    
    for node in fluid_model_part.Nodes:
	node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)   
	node.SetSolutionStepValue(DENSITY,0, density)   
    solver.Initialize()
    
elif(SolverType == "Quasi_Inc_Linear_Pressure"): 
    solver = ulf_fsi_inc.ULF_FSISolver(outputfile1, fluid_model_part, structure_model_part, combined_model_part, compute_reactions, box_corner1, box_corner2, domain_size, add_nodes, bulk_modulus, density)
    solver.alpha_shape = fluid_ulf_var.alpha_shape;
    solver.echo_level = 2;        
    for node in fluid_model_part.Nodes:
	node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)   
	node.SetSolutionStepValue(DENSITY,0, density)   
    solver.Initialize()
print "fluid solver created"

if (fluid_ulf_var.FSI==1):
    if (fluid_ulf_var.domain_size==2):
	fluid_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
    elif (fluid_ulf_var.domain_size==3):
	fluid_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
    else:
	raise "Domain size error. It should be 2D or 3D"

for node in fluid_model_part.Nodes:
  if(node.GetSolutionStepValue(DENSITY) == 0.0):
    print "node ",node.Id," has zero density!"
    raise 'node with zero density found'
  if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
    print "node ",node.Id," has zero viscosity!"
    raise 'node with zero VISCOSITY found'  
  if(FSI==1):
    is_fsi_interf+=node.GetSolutionStepValue(IS_INTERFACE)

if (SolverType == "Incompressible_Modified_FracStep" and FSI==1):
    if (is_fsi_interf==0):
	raise 'For running FSI using the Modified Frac Step Solver you must prescribe IS_INTERFACE flag at the surface/outer contour of your structure'

##############################################################################################                   
#   (_)(_)(_)(_)            (_)(_)(_)(_)(_)         (_) _     _ (_)   
#    (_)      (_)_          (_)                     (_)(_)   (_)(_)   
#    (_)        (_)         (_) _  _                (_) (_)_(_) (_)   
#    (_)        (_)         (_)(_)(_)               (_)   (_)   (_)   
#    (_)       _(_)         (_)                     (_)         (_)   
#    (_)_  _  (_)           (_) _  _  _  _          (_)         (_)   
#   (_)(_)(_)(_)            (_)(_)(_)(_)(_)         (_)         (_)   
################################################################################################   
#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
solid_model_part.SetBufferSize(2)
##adding dofs
SolverStrategy.AddDofs(solid_model_part)

#creating a DEM_solver object

DEM_solver = SolverStrategy.ExplicitStrategy(solid_model_part, domain_size) #here, solver variables initialize as default

if(ConcreteTestOption =="ON"):
  contact_model_part = DEM_solver.contact_model_part   
  
ProcGiDSolverTransfer(solid_model_part,DEM_solver)
  
DEM_solver.Initialize()
DEM_dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)
DEM_dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

if(ConcreteTestOption =="ON"):

  ProcModelData(solid_model_part,DEM_solver)       # calculates the mean number of neighbours the mean radius, etc..
  ProcListDefinition(solid_model_part,DEM_solver)  # defines the lists where we measure forces
  ProcSkinAndPressure(solid_model_part,DEM_solver)       # defines the skin and applies the pressure

if (CriticalTimeOption =="ON"):
  DEM_solver.Initial_Critical_Time() 
  if (DEM_dt!=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)):
    print("WARNING: Delta time has been modifyed to the critical one")
    DEM_dt=solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

#initializations
DEM_time = 0.0
DEM_step = 0
DEM_time_old_print = 0.0
initial_pr_time = timer.clock()
initial_real_time = timer.time()

print ('\n'+'Calculation starts at instant: ' + str(initial_pr_time)+'\n')

total_steps_expected = int(final_time/DEM_dt)
print ('Total number of TIME STEPs expected in the calculation is: ' + str(total_steps_expected) + ' if time step is kept ' +'\n' )

#paths:
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
results = open('results.txt','w') #file to export some results
summary_results = open('summary_results.txt','w')

forcelist = []
forcelist2 = []
timelist = []
displacementlist = []

os.chdir(list_path)
DEM_multifile = open(problem_name+'_all'+'.post.lst','w')
multifile_5 = open(problem_name+'_5'+'.post.lst','w')
multifile_10 = open(problem_name+'_10'+'.post.lst','w')
multifile_50 = open(problem_name+'_50'+'.post.lst','w')

DEM_multifile.write('Multiple\n')
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
contact_model_part = DEM_solver.contact_model_part   

os.chdir(post_path)
if(Multifile == "single_file"):
  DEM_gid_io.InitializeMesh(0.0)
  DEM_gid_io.WriteMesh(contact_model_part.GetMesh());
  DEM_gid_io.FinalizeMesh()
  DEM_gid_io.InitializeResults(0.0, contact_model_part.GetMesh()); 
  DEM_gid_io.InitializeMesh(0.0)
  DEM_gid_io.WriteSphereMesh(solid_model_part.GetMesh())
  DEM_gid_io.FinalizeMesh()
  DEM_gid_io.InitializeResults(0.0, solid_model_part.GetMesh()); 

os.chdir(main_path)
# for plotting the graph:     
velocity_node_y = 0.0
    
for node in force_measurement:
    velocity_node_y = node.GetSolutionStepValue(VELOCITY_Y,0) #Applied velocity during the uniaxial compression test

#done=False  #flag for the end of the confinement  
 
################################################################################################
#                #
#               ###                
#              #####              
#             #######
#            #########
#               ###
#               ###     
#################################################################################################

#settings to be changed
Dt = fluid_ulf_var.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01
final_time = fluid_ulf_var.max_time
output_step = fluid_ulf_var.output_step
safety_factor = 0.5 #you should put a safety factor ;-)!!!
next_output_time = output_step
time = 0.0
step = 0
inlet_vel = Vector(3)
if (fluid_ulf_var.lagrangian_nodes_inlet==1):    
    for node in fluid_model_part.Nodes:
	if (node.GetSolutionStepValue(IS_LAGRANGIAN_INLET)==1):
	    inlet_vel=node.GetSolutionStepValue(VELOCITY,0)
	    print "Lagrangian Inlet(s) Velocity  is ", inlet_vel
	    break
else:
    inlet_vel[0]=0.0
    inlet_vel[1]=0.0
    inlet_vel[2]=0.0

dummy=LagrangianInletProcess(fluid_model_part, 0.0, inlet_vel)
while (time < final_time):
    print "IN THE LOOP"
    step = step+1
    print "the new step is ",step
    print time
    if(step <= 3):
        new_Dt = 0.00000001;
        time = time + new_Dt*safety_factor

    #solving the fluid problem
    if(step > 3):
        print "the domain size is ", domain_size
        print "the Dt: ", Dt
        new_Dt = solver.EstimateDeltaTime(Dt, domain_size)
        print "the new Dt: ", new_Dt
        time = time + new_Dt*safety_factor

        combined_model_part.CloneTimeStep(time)        
        solver.Solve(dummy)
               
        print "after completing the solution"

        if(time > next_output_time):
    
            file_name = input_file_name
            file_name = file_name + str(time)

            gid_io.InitializeMesh( time );
            gid_io.WriteNodeMesh((combined_model_part).GetMesh());
            gid_io.WriteMesh((combined_model_part).GetMesh());
            gid_io.FinalizeMesh();

            gid_io.InitializeResults(time, (combined_model_part).GetMesh());
    
            gid_io.WriteNodalResults(DISPLACEMENT, combined_model_part.Nodes, time, 0);            
            gid_io.WriteNodalResults(VELOCITY, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(PRESSURE, (combined_model_part).Nodes, time, 0);
            if (compute_reactions==1):
		gid_io.WriteNodalResults(REACTION, (combined_model_part).Nodes, time, 0);
	    
            if (fluid_ulf_var.lagrangian_nodes_inlet==1):
                gid_io.WriteNodalResults(IS_LAGRANGIAN_INLET, (combined_model_part).Nodes, time, 0);                    
            
            gid_io.Flush()
            #gid_io.CloseResultFile();
            gid_io.FinalizeResults()
            next_output_time = next_output_time  + output_step;
##############################################################################################                   
#   (_)(_)(_)(_)            (_)(_)(_)(_)(_)         (_) _     _ (_)   
#    (_)      (_)_          (_)                     (_)(_)   (_)(_)   
#    (_)        (_)         (_) _  _                (_) (_)_(_) (_)   
#    (_)        (_)         (_)(_)(_)               (_)   (_)   (_)   
#    (_)       _(_)         (_)                     (_)         (_)   
#    (_)_  _  (_)           (_) _  _  _  _          (_)         (_)   
#   (_)(_)(_)(_)            (_)(_)(_)(_)(_)         (_)         (_)   
################################################################################################   
    DEM_dt = solid_model_part.ProcessInfo.GetValue(DELTA_TIME) #possible modifications of DELTA_TIME
    DEM_time += DEM_dt
    solid_model_part.CloneTimeStep(DEM_time)
    solid_model_part.ProcessInfo[TIME_STEPS] = DEM_step
        
    #printing forces in a file   
    total_force=0
    force_node= 0   
    os.chdir(data_and_results)   
    for node in force_measurement:
	force_node = node.GetSolutionStepValue(RHS,0)
	force_node_x = node.GetSolutionStepValue(RHS,0)[0]
	force_node_y = node.GetSolutionStepValue(RHS,0)[1]
	force_node_z = node.GetSolutionStepValue(RHS,0)[2]
	results.write(str(node.Id)+"  "+str(DEM_step)+"  "+str(force_node_y)+'\n')
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
    for node in inf_layer:
	
		force_node = node.GetSolutionStepValue(RHS,0)
		force_node_x = node.GetSolutionStepValue(RHS,0)[0]
		force_node_y = node.GetSolutionStepValue(RHS,0)[1]
		force_node_z = node.GetSolutionStepValue(RHS,0)[2]
		total_force += force_node_y
       
    #writing lists to be printed
    forcelist2.append(total_force)  
    summary_results.write(str(DEM_step)+"  "+str(total_force)+'\n')
    
    os.chdir(main_path)    
    #Dissable confinement   
    #if(ConcreteTestOption==True):
    
	  #if( (time > final_time*0.1) and (done==False)):
		#done=True;
		#for element in skin_list:
		  #element.SetValue(APPLIED_FORCE,(0,0,0))
		
		#print("Confinement finished at time "+str(time))
       
    DEM_solver.Solve() 
    incremental_time = (timer.time()-initial_real_time)- prev_time
    if (incremental_time > control_time):       
      percentage = 100.0*(float(DEM_step)/total_steps_expected)
      print 'Real time calculation: ' + str(timer.time()-initial_real_time) 
      print 'Percentage Completed: ' +str(percentage) + ' %' 
      print "TIME STEP = " + str(DEM_step) + '\n'    	
      prev_time = (timer.time()-initial_real_time)
    
    if ( (timer.time()-initial_real_time > 60.0) and cond==0):	
      cond=1	
      estimation_time=60.0*(total_steps_expected/DEM_step) #seconds	
      print('the total calculation estimated time is '+str(estimation_time)+'seconds.'+'\n')
      print('in minutes :'+str(estimation_time/60)+'min.'+'\n')
      print('in hours :'+str((estimation_time/60)/60)+'hrs.'+'\n')
      print('in days :'+str(((estimation_time/60)/60)/24)+'days.'+'\n')		
      if (((estimation_time/60)/60)/24 > 2.0):
	print('WARNING!!!:       VERY LASTING CALCULATION'+'\n')
	   	      
    os.chdir(list_path)    
    DEM_multifile.write(problem_name+'_'+str(time)+'.post.bin\n')   
    os.chdir(main_path)

  ##############     GiD IO        ################################################################################
    
    time_to_print = DEM_time - DEM_time_old_print
    #print str(time)    
    if(time_to_print >= output_dt):   
	os.chdir(graphs_path)
	#Drawing graph stress_strain:
	if(ConcreteTestOption =="ON"):
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
	  DEM_gid_io.InitializeMesh(time)
	  DEM_gid_io.WriteSphereMesh(solid_model_part.GetMesh())
	  DEM_gid_io.FinalizeMesh()
	  DEM_gid_io.InitializeResults(time, solid_model_part.GetMesh()); 
	  DEM_gid_io.InitializeMesh(time)
	  DEM_gid_io.WriteMesh(contact_model_part.GetMesh());
	  DEM_gid_io.FinalizeMesh()
	  DEM_gid_io.InitializeResults(time, contact_model_part.GetMesh()); 
 
	##########PRINTING VARIABLES############
	
	ProcPrintingVariables(DEM_gid_io,solid_model_part,contact_model_part,time)  
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
	  DEM_gid_io.FinalizeResults()
	
	os.chdir(main_path)                  
	DEM_time_old_print = time
   
    #End of print loop

    os.chdir(main_path)    
    graph_export.write(str(strain)+"  "+str(total_stress)+'\n')    
    DEM_step += 1
 
################################################################################################
#                #
#               ###                
#              #####              
#             #######
#            #########
#               ###
#               ###     
#################################################################################################
 

 
##############################################################################################                   
#   (_)(_)(_)(_)            (_)(_)(_)(_)(_)         (_) _     _ (_)   
#    (_)      (_)_          (_)                     (_)(_)   (_)(_)   
#    (_)        (_)         (_) _  _                (_) (_)_(_) (_)   
#    (_)        (_)         (_)(_)(_)               (_)   (_)   (_)   
#    (_)       _(_)         (_)                     (_)         (_)   
#    (_)_  _  (_)           (_) _  _  _  _          (_)         (_)   
#   (_)(_)(_)(_)            (_)(_)(_)(_)(_)         (_)         (_)   
################################################################################################  
if(Multifile == "single_file"):
  DEM_gid_io.FinalizeResults()

os.chdir(data_and_results)
graph_export.close() 
results.close()
summary_results.close()

os.chdir(list_path)
DEM_multifile.close()
multifile_5.close()
multifile_10.close()
multifile_50.close()

#plots
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
print "ANALYSIS COMPLETED" 
################################################################################################
#                #
#               ###                
#              #####              
#             #######
#            #########
#               ###
#               ###     
#################################################################################################