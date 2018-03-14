from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ChimeraApplication import *
import csv 

######################################################################################
######################################################################################
######################################################################################
##PARSING THE PARAMETERS
#import define_output

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

parameter_file = open("ProjectParameters_backup.json",'r')
ProjectParameters = Parameters( parameter_file.read())


## Fluid model part definition 
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

## Solver construction
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])
################################################DISTANCE VARIABLE######################################
#main_model_part.AddNodalSolutionStepVariable(FUEL)
main_model_part.AddNodalSolutionStepVariable(DISTANCE)
main_model_part.AddNodalSolutionStepVariable(NORMAL)
main_model_part.AddNodalSolutionStepVariable(NODAL_MASS)
solver.AddVariables()

## Read the model - note that SetBufferSize is done here
solver.ImportModelPart()

## Add AddDofs
solver.AddDofs()

## Modelparts
background = main_model_part.GetSubModelPart("GENERIC_background")
patch = main_model_part.GetSubModelPart("GENERIC_patch")
patchBoundary = main_model_part.GetSubModelPart("GENERIC_patchBoundary")
#background.AddNodalSolutionStepVariable(DISTANCE)
#patch.AddNodalSolutionStepVariable(DISTANCE)

## Initialize GiD  I/O
from gid_output_process import GiDOutputProcess
gid_output_backGround = GiDOutputProcess(background,
                              "BackgroundSolution" ,
                              ProjectParameters["output_configuration"])

gid_output_patch = GiDOutputProcess(patch,
                              "PatchSolution" ,
                              ProjectParameters["output_configuration"])

gid_output_backGround.ExecuteInitialize()
gid_output_patch.ExecuteInitialize()

##here all of the allocation of the strategies etc is done
solver.Initialize()


##TODO: replace MODEL for the Kratos one ASAP
## Get the list of the skin submodel parts in the object Model
for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
    skin_part_name = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
    Model.update({skin_part_name: main_model_part.GetSubModelPart(skin_part_name)})

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
# Note that the conditions are firstly constructed. Otherwise, they may overwrite the BCs information.
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["initial_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["gravity"] )

## Processes initialization
for process in list_of_processes:
    process.ExecuteInitialize()


# Initialise VtkOutput
vtkOutput = VtkOutput(main_model_part,"nnn",ProjectParameters["output_configuration"])

#TODO: think if there is a better way to do this
#fluid_model_part = solver.GetComputeModelPart()

## Stepping and time settings
Dt = ProjectParameters["problem_data"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

time = 0.0
step = 0
out = 0.0

gid_output_backGround.ExecuteBeforeSolutionLoop()
gid_output_patch.ExecuteBeforeSolutionLoop()

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
    
    
file_name_patch = 'patch.csv'  
with open (file_name_patch,'w') as file:
  file.write('Time,massflux,PRESSURE,VELOCITY_X,VELOCITY_Y\n') 
  

file_name_background = 'background.csv'  
with open (file_name_background,'w') as file:
  file.write('Time,massflux,PRESSURE,VELOCITY_X,VELOCITY_Y\n') 
    

def CalculateMassFlux(model_part,time,filename,node_r):
  flux =0
  
  for node in model_part.Nodes:
    vel_x = node.GetSolutionStepValue(VELOCITY_X,0)
    vel_y = node.GetSolutionStepValue(VELOCITY_Y,0)
    n_x = node.GetSolutionStepValue(NORMAL,0)[0]
    n_y = node.GetSolutionStepValue(NORMAL,0)[1]
    flux = flux + vel_x*n_x+vel_y*n_y 
    pressure = model_part.GetNode(node_r).GetSolutionStepValue(PRESSURE,0)
    velocity_x = model_part.GetNode(node_r).GetSolutionStepValue(VELOCITY_X,0)
    velocity_y = model_part.GetNode(node_r).GetSolutionStepValue(VELOCITY_Y,0)
    
  with open(filename, 'a') as file:
    file.write(str(time)+","+str(flux)+","+str(pressure)+','+str(velocity_x)+","+str(velocity_y)+'\n')
    
node_patch = 1264
node_background = 428
    

ChimeraImplementer = CustomApplyChimeraUsingMpcProcess2d(main_model_part,background,patch,59)
ChimeraImplementer.ApplyChimeraUsingMpc2d(patchBoundary,"NearestElement")


while(time <= end_time):

    time = time + Dt
    step = step + 1 
    main_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)
    
    if(step >= 3):
        #print("Nr of processes",len(list_of_processes))
        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()
            #print("process ", process)
        
            
        gid_output_backGround.ExecuteInitializeSolutionStep()
        gid_output_patch.ExecuteInitializeSolutionStep()
        
        
        
        
        solver.Solve()
        #movePatch(patch)
        #print("exiting the process loop")
        
        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()
            
        gid_output_backGround.ExecuteFinalizeSolutionStep()
        gid_output_patch.ExecuteFinalizeSolutionStep()

        #TODO: decide if it shall be done only when output is processed or not
        for process in list_of_processes:
            process.ExecuteBeforeOutputStep()
    
        if gid_output_backGround.IsOutputStep():
            gid_output_backGround.PrintOutput()
        if gid_output_patch.IsOutputStep():
            gid_output_patch.PrintOutput()
            vtkOutput.PrintOutput()        
        for process in list_of_processes:
            process.ExecuteAfterOutputStep()
            
        CalculateMassFlux(patch,time,file_name_patch,node_patch)
        CalculateMassFlux(background,time,file_name_background,node_background)
      
        

        out = out + Dt

for process in list_of_processes:
    process.ExecuteFinalize()
    
gid_output_backGround.ExecuteFinalize()
gid_output_patch.ExecuteFinalize()














