from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.FluidRveLagrangeMultipliersApplication import *
from math import sqrt
from math import sin
######################################################################################
######################################################################################
######################################################################################
##PARSING THE PARAMETERS
#import define_output

parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()

## Import KratosMPI if needed
if (parallel_type == "MPI"):
    import KratosMultiphysics.mpi as KratosMPI

## Fluid model part definition
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

## Solver construction
import navier_stokes_solver_vmsmonolithic_modified
solver = navier_stokes_solver_vmsmonolithic_modified.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

solver.AddVariables()
#NEW LINES: adding normal and nodal area:
main_model_part.AddNodalSolutionStepVariable(NORMAL)
main_model_part.AddNodalSolutionStepVariable(NODAL_AREA)
main_model_part.AddNodalSolutionStepVariable(LAGRANGE_MULTIPLIER_VELOCITY)
main_model_part.AddNodalSolutionStepVariable(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS)
main_model_part.AddNodalSolutionStepVariable(NODE_PAIR_X_COMPONENT)
main_model_part.AddNodalSolutionStepVariable(NODE_PAIR_Y_COMPONENT)
main_model_part.AddNodalSolutionStepVariable(NODE_PAIR_PRESSURE)
main_model_part.AddNodalSolutionStepVariable(NODE_PAIR_X_COMPONENT_ANTIPERIODIC)
main_model_part.AddNodalSolutionStepVariable(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC)
main_model_part.AddNodalSolutionStepVariable(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY)
main_model_part.AddNodalSolutionStepVariable(REYNOLDS_STRESS_2D)## Read the model - note that SetBufferSize is done here
solver.ImportModelPart()

## Add AddDofs
solver.AddDofs()



#NEW LINES: adding new dofs
node_1=True
node_2=True
node_i=-1

for node in main_model_part.Nodes:
   node.AddDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY)
for node in main_model_part.Nodes:
   node_i=node_i+1
   if node_i<2:
     node.AddDof(LAGRANGE_MULTIPLIER_VELOCITY_X)
     node.AddDof(LAGRANGE_MULTIPLIER_VELOCITY_Y)
     if (main_model_part.ProcessInfo.GetValue(DOMAIN_SIZE)==3):
             node.AddDof(LAGRANGE_MULTIPLIER_VELOCITY_Z)
     node.AddDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X)
     node.AddDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y)
     node.AddDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Z)
   else:
      break;


##NEW LINES: deleting condittions
for condition in main_model_part.Conditions:
    condition.Set(TO_ERASE,True)
condition_erase_process = ConditionEraseProcess(main_model_part)
condition_erase_process.Execute()
##NEW LINES: removing imposed velocities
for node in main_model_part.Nodes:
    if node.IsFixed(VELOCITY_X):
          node.Free(VELOCITY_X)
    if node.IsFixed(VELOCITY_Y):
          node.Free(VELOCITY_Y)
    if node.IsFixed(VELOCITY_Z):
          node.Free(VELOCITY_Z)
    if node.IsFixed(PRESSURE):
          node.Free(PRESSURE)
    '''
    if node.Y<0.001 or node.Y>0.999:
          node.Fix(VELOCITY_X)
          node.Fix(VELOCITY_Y)
          node.SetSolutionStepValue(VELOCITY_X,1.0)
    if node.X<0.001:
          node.Fix(VELOCITY_Y)
          node.Fix(VELOCITY_X)
          node.SetSolutionStepValue(VELOCITY_X,1.0)
    if node.X>0.999:
          node.Fix(VELOCITY_Y)
 
    if node.X>0.999:
          node.Fix(PRESSURE)
    '''
    #if node.Y<0.0001:
    #      node.Fix(VELOCITY_X)
    #node.Fix(VELOCITY_Y)

viscosity = 0.01
density = 1.0*1.0;

first_node=True
node_i=-1
for node in main_model_part.Nodes:
      node_i=node_i+1
      node.SetSolutionStepValue(VISCOSITY,viscosity)
      node.SetSolutionStepValue(DENSITY,density)
      node.SetSolutionStepValue(VELOCITY_X,0.0)
      y_coord=node.Y
      node.Free(VELOCITY_X)
      node.Free(VELOCITY_Y)
      '''
      if node.Id==111:
          node.Fix(PRESSURE)
          node.SetSolutionStepValue(PRESSURE,0.0)
      if node.Id==347:
          node.Fix(PRESSURE)
          node.SetSolutionStepValue(PRESSURE,1.0)
      '''
      #if node.X<0.001:
      #   node.Fix(VELOCITY_X)
      #if node.Y<0.001:
      #   node.Fix(VELOCITY_Y)
      #node.SetSolutionStepValue(VELOCITY_X,(node.Y-0.5))
      #node.SetSolutionStepValue(VELOCITY_Y,-(node.X-0.5))
      #node.SetSolutionStepValue(VELOCITY_X,1/1.41)
      #node.SetSolutionStepValue(VELOCITY_Y,1.0)
      #node.SetSolutionStepValue(VELOCITY_X,0.0+(node.Y-0.5)*10.0)
      #node.SetSolutionStepValue(VELOCITY_Y,0.0+(node.X-0.5)*10.0)
      #node.SetSolutionStepValue(VELOCITY_Y,0.25-(node.X-0.5)*(node.X-0.5))
      #node.SetSolutionStepValue(VELOCITY_X,0.25-(node.Y-0.5)*(node.Y-0.5))
      if first_node:
          node.Fix(PRESSURE)
          first_node=False
      if node_i<2:
          1
          #node.Fix(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y)
          node.Fix(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Z)
          #node.Fix(LAGRANGE_MULTIPLIER_VELOCITY_X)
          #node.Fix(LAGRANGE_MULTIPLIER_VELOCITY_Y)
#NEW LINES: finding normal and nodal area-
normal_tools = BodyNormalCalculationUtils()
normal_tools.CalculateBodyNormals(main_model_part,main_model_part.ProcessInfo.GetValue(DOMAIN_SIZE))  

#NEW LINES: setting the target velocity:
#NEW LINES: setting the target velocity gradient:
node_i=-1
for node in main_model_part.Nodes:
  node_i=node_i+1
  if node_i==0:
      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_X,0.0)#1/1.41)
      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_Y,0.0)#1/1.41)
      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_Z,0.0)
      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X,160.00)
      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y,0.0)
      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Z,0.0)
  elif node_i==1:
      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_X,0.0)#1/1.41)
      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_Y,0.0)#1/1.41)
      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_Z,0.0)
      #node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X,-10.0)
      #node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y,-10.0)
      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Z,0.0)
  else:
    break;
#NEW LINES: defining the imposed strains:
strain=Vector(3)
strain[0]=0.0
strain[1]=0.0
strain[2]=0.0
main_model_part.ProcessInfo.SetValue(STRAIN,strain)


#NEW LINES: adding the lagrange multipliers
add_lagrange_multipliers_util = AddMeanVelocityLagrangeMultiplierConditions2D(main_model_part)
add_lagrange_multipliers_util.AddThem(-1000.0,1000.0,-1000.0,1000.0)
add_lagrange_multipliers_util_grad = AddVelocityGradientsLagrangeMultiplierConditions2D(main_model_part)
####add_lagrange_multipliers_util_grad.AddThem(-1000.0,1000.0,-1000.0,1000.0)
y_division=-0.5
#add_lagrange_multipliers_util_grad.AddThem2(-1000.0,1000.0,-1000.0,1000.0,y_division)
add_periodic_conditions = AddPeriodicConditionsNormalOnly2D(main_model_part)
#add_periodic_conditions.AddThem(0.0001, 0.999 , 0.0001 ,0.999) 
add_periodic_conditions.AddThemWithTangentInversePeriodicity(0.0001, 0.999 , 0.0001 ,0.999) 


vel = sqrt ( main_model_part.ProcessInfo.GetValue(LAGRANGE_MULTIPLIER_VELOCITY_X) * main_model_part.ProcessInfo.GetValue(LAGRANGE_MULTIPLIER_VELOCITY_X) + main_model_part.ProcessInfo.GetValue(LAGRANGE_MULTIPLIER_VELOCITY_Y)*main_model_part.ProcessInfo.GetValue(LAGRANGE_MULTIPLIER_VELOCITY_Y))
lenght = 1.0
rho = 1.0

#grad_total = main_model_part.ProcessInfo.GetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X) +  main_model_part.ProcessInfo.GetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y)

print("el reynolds es:", vel * rho *lenght /viscosity)
#print("el Idelsohn es:", grad_total * rho * lenght*lenght /viscosity)
#print(hola)


##here all of the allocation of the strategies etc is done
solver.Initialize()




##TODO: replace MODEL for the Kratos one ASAP
## Get the list of the skin submodel parts in the object Model
for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
    skin_part_name = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
    Model.update({skin_part_name: main_model_part.GetSubModelPart(skin_part_name)})

#~ ## Get the list of the no-skin submodel parts in the object Model
#~ for i in range(ProjectParameters["solver_settings"]["no_skin_parts"].size()):
    #~ no_skin_part_name = ProjectParameters["solver_settings"]["no_skin_parts"][i].GetString()
    #~ Model.update({no_skin_part_name: main_model_part.GetSubModelPart(no_skin_part_name)})

## Get the list of the initial conditions submodel parts in the object Model
for i in range(ProjectParameters["initial_conditions_process_list"].size()):
    initial_cond_part_name = ProjectParameters["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
    Model.update({initial_cond_part_name: main_model_part.GetSubModelPart(initial_cond_part_name)})

## Get the gravity submodel part in the object Model
for i in range(ProjectParameters["gravity"].size()):
    gravity_part_name = ProjectParameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
    Model.update({gravity_part_name: main_model_part.GetSubModelPart(gravity_part_name)})


## Processes construction
import process_factory
# "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity)
# Note 1: gravity is firstly constructed. Outlet process might need its information.
# Note 2: conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["gravity"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["initial_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )

## Processes initialization
for process in list_of_processes:
    process.ExecuteInitialize()

#TODO: think if there is a better way to do this
fluid_model_part = solver.GetComputingModelPart()



gid_mode = GiDPostMode.GiD_PostBinary#PostAscii
multifile = MultiFileFlag.SingleFile #MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
"""Result file"""
gid_io = GidIO("caso_email_2905_visc0625_den1_grad40_with_reytensor_noanti",gid_mode,multifile,deformed_mesh_flag,write_conditions)
## Initialize GiD  I/O
gid_io.InitializeMesh( 0.0 );
gid_io.WriteMesh((fluid_model_part).GetMesh());
gid_io.FinalizeMesh()
"""Result initialization"""
gid_io.InitializeResults(0.0,(fluid_model_part).GetMesh())





for node in main_model_part.Nodes:
   1
   #if node.Y>0.5:
   #   node.SetSolutionStepValue(VELOCITY_X,(node.Y-0.75)*10.0)
   #else:
   #    node.SetSolutionStepValue(VELOCITY_X,-(node.Y-0.25)*10.0)
   
   node.SetSolutionStepValue(VELOCITY_X,(node.Y-0.5)*80.0*0.25)
   node.SetSolutionStepValue(VELOCITY_Y,(node.X-0.5)*80.0*0.25)
   #node.SetSolutionStepValue(VELOCITY_X,0.01*(sin(node.Y*2.0*3.14159)*40.0-(node.Y-0.5)*80.0))
   #node.SetSolutionStepValue(VELOCITY_Y,0.01*(sin(node.X*2.0*3.14159)*40.0-(node.X-0.5)*80.0))
   #node.SetSolutionStepValue(VELOCITY_Y,(node.X-0.5)*10.0)
   #if node.X>0.4 and node.X<0.6: # and node.Y>0.3 and node.Y<0.7:
   #   node.SetSolutionStepValue(VELOCITY_Y,4.0)
   #node.SetSolutionStepValue(VELOCITY_X,1,(node.Y-0.5)*5.0)
   #node.SetSolutionStepValue(VELOCITY_Y,(node.X-0.5)*5.0)
   #node.SetSolutionStepValue(VELOCITY_Y,1,(node.X-0.5)*5.0)
   #node.SetSolutionStepValue(BODY_FORCE_X,1.001)
   #node.SetSolutionStepValue(VELOCITY_X,5.0)
   #if node.X>0.6 and node.X<0.8 and node.Y>0.4 and node.Y<0.6:
   #         node.SetSolutionStepValue(VELOCITY_Y,1.0)

## Stepping and time settings
Dt = ProjectParameters["problem_data"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

time = 0.0
step = 0
out = 0.0


for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

first_steps_dt= Dt*1
end_time = end_time + first_steps_dt*3.0

gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,0.0,0)

out_step=0


#NEW LINES: postprocess stress:
stress_util = ShearTermsComputationUtility2D(fluid_model_part)
first_row_stress = Array3()
second_row_stress = Array3()
stress_util.ComputeTau(first_row_stress,second_row_stress)
out_file1 = open("stresses_0625_half_visc.dat",'w')
out_line1 = '#time	#TAUxx	#TAUxy	#TAUyx	#TAUyy'+'\n'
out_file1.write(out_line1)
out_line1 = '0.0'+'	'+str(first_row_stress[0])+'	'+str(first_row_stress[1])+'	'+str(second_row_stress[0])+'	'+str(second_row_stress[1])+'\n'
out_file1.write(out_line1)
out_file1.flush()

'''
Dt=1.0
time=1.0
if True:
        for node in main_model_part.Nodes:
           vel = node.GetSolutionStepValue(VELOCITY)
           pos = Vector(3)
           pos[0] = node.Y-0.5
           pos[1] = node.X-0.5
           pos[2] = 0.0
           vel_perturbation = vel  - 80.0 * pos
           stress_0 = vel_perturbation[0]*vel_perturbation[0]
           stress_1 = vel_perturbation[1]*vel_perturbation[1]
           stress_2 = vel_perturbation[0]*vel_perturbation[1]
           old_stresses = node.GetSolutionStepValue(REYNOLDS_STRESS_2D)
           new_stress_0 = ( old_stresses[0] * (time-Dt) + stress_0*Dt ) / time
           new_stress_1 = ( old_stresses[1] * (time-Dt) + stress_1*Dt ) / time
           new_stress_2 = ( old_stresses[2] * (time-Dt) + stress_2*Dt ) / time
           node.SetSolutionStepValue(REYNOLDS_STRESS_2D_X,new_stress_0)
           node.SetSolutionStepValue(REYNOLDS_STRESS_2D_Y,new_stress_1)
           node.SetSolutionStepValue(REYNOLDS_STRESS_2D_Z,new_stress_2)
           print("vel=",vel)
           print("pos=",pos)
           print("vel_per=",vel_perturbation)
        print(asfasd)
'''
while(time <= end_time):
    out_step=out_step+1
    if step<=2:
      time = time + first_steps_dt
    else:
       time = time + Dt
    step = step + 1
    main_model_part.CloneTimeStep(time)

    if (parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
        print("STEP = ", step)
        print("TIME = ", time)

    if(step >= 3):
        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()
        '''
        if time<2.0:
            main_model_part.ProcessInfo.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X,time*0.5*10.0)
            main_model_part.ProcessInfo.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y,-time*0.5*10.0)
        else:
            main_model_part.ProcessInfo.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X,10.0)
            main_model_part.ProcessInfo.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y,0.0)
            main_model_part.ProcessInfo.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Z,0.0)
        '''
        solver.Solve()

        #if time>8.0:
        #   for node in main_model_part.Nodes:
        #      node.SetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X,140.0)
        #      break

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()


        #saving the reynolds stress.
        for node in main_model_part.Nodes:
           vel = node.GetSolutionStepValue(VELOCITY)
           pos = Vector(3)
           pos[0] = node.Y-0.5
           pos[1] = node.X-0.5
           pos[2] = 0.0
           vel_perturbation = vel  - 20.0 * pos
           stress_0 = vel_perturbation[0]*vel_perturbation[0]
           stress_1 = vel_perturbation[1]*vel_perturbation[1]
           stress_2 = vel_perturbation[0]*vel_perturbation[1]
           old_stresses = node.GetSolutionStepValue(REYNOLDS_STRESS_2D)
           new_stress_0 = ( old_stresses[0] * (time-Dt) + stress_0*Dt ) / time
           new_stress_1 = ( old_stresses[1] * (time-Dt) + stress_1*Dt ) / time
           new_stress_2 = ( old_stresses[2] * (time-Dt) + stress_2*Dt ) / time
           node.SetSolutionStepValue(REYNOLDS_STRESS_2D_X,new_stress_0)
           node.SetSolutionStepValue(REYNOLDS_STRESS_2D_Y,new_stress_1)
           node.SetSolutionStepValue(REYNOLDS_STRESS_2D_Z,new_stress_2)
           #print(vel_perturbation)

        #if gid_output.IsOutputStep():
        #    gid_output.PrintOutput()
        if out_step>=50:
         stress_util.ComputeTau(first_row_stress,second_row_stress)
         out_line1 = str(time)+'	'+str(first_row_stress[0])+'	'+str(first_row_stress[1])+'	'+str(second_row_stress[0])+'	'+str(second_row_stress[1])+'\n'
         out_file1.write(out_line1)
         out_file1.flush()
         gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
         gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
         gid_io.WriteNodalResults(LAGRANGE_MULTIPLIER_VELOCITY,fluid_model_part.Nodes,time,0)
         gid_io.WriteNodalResults(REYNOLDS_STRESS_2D,fluid_model_part.Nodes,time,0)
         #gid_io.WriteNodalResults(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS,fluid_model_part.Nodes,time,0)
         #gid_io.WriteNodalResults(NODE_PAIR_X_COMPONENT,fluid_model_part.Nodes,time,0)
         #gid_io.WriteNodalResults(NODE_PAIR_Y_COMPONENT,fluid_model_part.Nodes,time,0)
         #gid_io.WriteNodalResults(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY,fluid_model_part.Nodes,time,0)
         #gid_io.WriteNodalResults(VISCOSITY,fluid_model_part.Nodes,time,0)
         print("pritinting")
         out_step=0        
 
        #saving the reynolds stress.


        for process in list_of_processes:
            process.ExecuteAfterOutputStep()

        out = out + Dt

for process in list_of_processes:
    process.ExecuteFinalize()

gid_io.FinalizeResults()
