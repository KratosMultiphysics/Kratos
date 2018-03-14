from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
# ADI BEGIN
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ChimeraApplication import *
import math


patchNameParameterFileMap = {"square":"ProjectParameters_innerMesh.json"}
backGroundParameterFile = "ProjectParameters.json"
patchNameParameterMap = {}
patchNameModelPartMap = {}
patchNameModelMap = {}
patchNameSolverMap = {}
patchNameGidOutPutMap = {}
patchNameListProcessesMap = {}

parameter_file = open(backGroundParameterFile,'r')
backGroundProjectParameters = Parameters( parameter_file.read())

for patchName, patchParameterFile in patchNameParameterFileMap.items():
    parameter_file = open(patchParameterFile,'r')
    ProjectParameters = Parameters( parameter_file.read())
    patchNameParameterMap[patchName] = ProjectParameters

## Fluid model part definition 
backGround_model_part = ModelPart(backGroundProjectParameters["problem_data"]["model_part_name"].GetString())
backGround_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, backGroundProjectParameters["problem_data"]["domain_size"].GetInt())

for patchName, patchProjectParameters in patchNameParameterMap.items():
    patchModelPart = ModelPart(patchProjectParameters["problem_data"]["model_part_name"].GetString())
    patchModelPart.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())
    patchNameModelPartMap[patchName] = patchModelPart
    
   
###TODO replace this "model" for real one once available
backGroundModel = {backGroundProjectParameters["problem_data"]["model_part_name"].GetString() : backGround_model_part}

for patchName, patchProjectParameters in patchNameParameterMap.items():
    patchModelPart = patchNameModelPartMap[patchName]
    patchModel = {patchProjectParameters["problem_data"]["model_part_name"].GetString() : patchModelPart}
    patchNameModelMap[patchName] = patchModel
    
print('################ :: Setting up background solver')
## Solver construction
backGround_solver_module = __import__(backGroundProjectParameters["solver_settings"]["solver_type"].GetString())
back_ground_solver = backGround_solver_module.CreateSolver(backGround_model_part, backGroundProjectParameters["solver_settings"])
back_ground_solver.AddVariables()
## Read the model - note that SetBufferSize is done here
back_ground_solver.ImportModelPart()
## Add AddDofs
back_ground_solver.AddDofs()


for patchName, patchProjectParameters in patchNameParameterMap.items():
    print('################ :: Setting up patch solver for :: ', patchName)
    patchModelPart = patchNameModelPartMap[patchName]    
    patch_solver_module = __import__(patchProjectParameters["solver_settings"]["solver_type"].GetString())
    patch_solver = patch_solver_module.CreateSolver(patchModelPart, patchProjectParameters["solver_settings"])
    patch_solver.AddVariables()
    ## Read the model - note that SetBufferSize is done here
    patch_solver.ImportModelPart()
    ## Add AddDofs
    patch_solver.AddDofs()
    patchNameSolverMap[patchName] = patch_solver


### Here we make a Chimera solver and add the back ground and chimera patch
### 
chimera_solver_module = __import__("partitioned_chimera_solver")
chimeraSolver = chimera_solver_module.CreateSolver()
chimeraSolver.AddBackgroundModelPart(backGround_model_part, back_ground_solver)
chimeraSolver.AddChimeraPatchModelParts(patchNameModelPartMap, patchNameSolverMap)
chimeraSolver.AddVariables()
#chimeraSolver.Initialize('Inlet3D_interface')
#print('###################### ')

## Initialize GiD  I/O
from gid_output_process import GiDOutputProcess
background_gid_output = GiDOutputProcess(back_ground_solver.GetComputeModelPart(),
                              backGroundProjectParameters["problem_data"]["problem_name"].GetString() ,
                              backGroundProjectParameters["output_configuration"])
background_gid_output.ExecuteInitialize()
##here all of the allocation of the strategies etc is done
back_ground_solver.Initialize()


for patchName, patchProjectParameters in patchNameParameterMap.items():
    patch_solver = patchNameSolverMap[patchName]    
    patch_gid_output = GiDOutputProcess(patch_solver.GetComputeModelPart(),
                              patchProjectParameters["problem_data"]["problem_name"].GetString() ,
                              patchProjectParameters["output_configuration"])
    patch_gid_output.ExecuteInitialize()
    patch_solver.Initialize()
    patchNameGidOutPutMap[patchName] = patch_gid_output


##TODO: replace MODEL for the Kratos one ASAP
## Get the list of the skin submodel parts in the object Model
for i in range(backGroundProjectParameters["solver_settings"]["skin_parts"].size()):
    skin_part_name = backGroundProjectParameters["solver_settings"]["skin_parts"][i].GetString()
    backGroundModel.update({skin_part_name: backGround_model_part.GetSubModelPart(skin_part_name)})

## Get the list of the initial conditions submodel parts in the object Model
for i in range(backGroundProjectParameters["initial_conditions_process_list"].size()):
    initial_cond_part_name = backGroundProjectParameters["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
    backGroundModel.update({initial_cond_part_name: backGround_model_part.GetSubModelPart(initial_cond_part_name)})
    
## Get the gravity submodel part in the object Model
for i in range(backGroundProjectParameters["gravity"].size()):   
    gravity_part_name = backGroundProjectParameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
    backGroundModel.update({gravity_part_name: backGround_model_part.GetSubModelPart(gravity_part_name)})    



for patchName, patchProjectParameters in patchNameParameterMap.items():
    patchModelPart = patchNameModelPartMap[patchName]    
    patchModel = patchNameModelMap[patchName]
    
    for i in range(patchProjectParameters["solver_settings"]["skin_parts"].size()):
        skin_part_name = patchProjectParameters["solver_settings"]["skin_parts"][i].GetString()
        patchModel.update({skin_part_name: patchModelPart.GetSubModelPart(skin_part_name)})

    ## Get the list of the initial conditions submodel parts in the object Model
    for i in range(patchProjectParameters["initial_conditions_process_list"].size()):
        initial_cond_part_name = patchProjectParameters["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
        patchModel.update({initial_cond_part_name: patchModelPart.GetSubModelPart(initial_cond_part_name)})
        
    ## Get the gravity submodel part in the object Model
    for i in range(patchProjectParameters["gravity"].size()):   
        gravity_part_name = patchProjectParameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
        patchModel.update({gravity_part_name: patchModelPart.GetSubModelPart(gravity_part_name)})
        


## Processes construction    
import process_factory
# "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity) 
# Note that the conditions are firstly constructed. Otherwise, they may overwrite the BCs information.
background_list_of_processes = process_factory.KratosProcessFactory(backGroundModel).ConstructListOfProcesses( backGroundProjectParameters["initial_conditions_process_list"] )
background_list_of_processes += process_factory.KratosProcessFactory(backGroundModel).ConstructListOfProcesses( backGroundProjectParameters["boundary_conditions_process_list"] )
background_list_of_processes += process_factory.KratosProcessFactory(backGroundModel).ConstructListOfProcesses( backGroundProjectParameters["gravity"] )
## Processes initialization
for process in background_list_of_processes:
    process.ExecuteInitialize()

for patchName, patchProjectParameters in patchNameParameterMap.items():
    patchModel = patchNameModelMap[patchName]    
    patch_list_of_processes = process_factory.KratosProcessFactory(patchModel).ConstructListOfProcesses( patchProjectParameters["initial_conditions_process_list"] )
    patch_list_of_processes += process_factory.KratosProcessFactory(patchModel).ConstructListOfProcesses( patchProjectParameters["boundary_conditions_process_list"] )
    patch_list_of_processes += process_factory.KratosProcessFactory(patchModel).ConstructListOfProcesses( patchProjectParameters["gravity"] )
    ## Processes initialization
    for process in patch_list_of_processes:
        process.ExecuteInitialize()
    patchNameListProcessesMap[patchName] = patch_list_of_processes
    
    
## Stepping and time settings
Dt = backGroundProjectParameters["problem_data"]["time_step"].GetDouble()
end_time = backGroundProjectParameters["problem_data"]["end_time"].GetDouble()

time = 0.0
step = 0
out = 0.0

background_gid_output.ExecuteBeforeSolutionLoop()
for process in background_list_of_processes:
    process.ExecuteBeforeSolutionLoop()
    
for patchName, patch_gid_output in patchNameGidOutPutMap.items():
    patch_gid_output.ExecuteBeforeSolutionLoop()
    patch_process_list = patchNameListProcessesMap[patchName]
    for process in patch_process_list:
        process.ExecuteBeforeSolutionLoop()    

backGround_model_part = patchNameModelPartMap['square']
background_list_of_processes = patchNameListProcessesMap['square']
background_gid_output = patchNameGidOutPutMap['square']
back_ground_solver = patchNameSolverMap['square']


while(time <= end_time):

    time = time + Dt
    step = step + 1
    backGround_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 3):
        for process in background_list_of_processes:
            process.ExecuteInitializeSolutionStep()
        
        background_gid_output.ExecuteInitializeSolutionStep()
        
        back_ground_solver.Solve()
        
        for process in background_list_of_processes:
            process.ExecuteFinalizeSolutionStep()
            
        background_gid_output.ExecuteFinalizeSolutionStep()

        #TODO: decide if it shall be done only when output is processed or not
        for process in background_list_of_processes:
            process.ExecuteBeforeOutputStep()
    
        if background_gid_output.IsOutputStep():
            background_gid_output.PrintOutput()
        
        for process in background_list_of_processes:
            process.ExecuteAfterOutputStep()

        out = out + Dt

for process in background_list_of_processes:
    process.ExecuteFinalize()
    
background_gid_output.ExecuteFinalize()
