from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
import KratosMultiphysics
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
# ADI BEGIN
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ChimeraApplication import *
import math
import numpy as np


patchNameParameterFileMap = {"square":"ProjectParameters_innerMesh.json"}
backGroundParameterFile = "ProjectParameters.json"
patchNameParameterMap = {}
patchNameModelPartMap = {}
patchNameModelMap = {}
patchNameSolverMap = {}
patchNameGidOutPutMap = {}
patchNameListProcessesMap = {}

''' ##########################################  '''
''' SETUP for CHIMERA ''' 
''' ##########################################  '''
### Here we make a Chimera solver and add the back ground and chimera patch
### 
chimera_solver_module = __import__("partitioned_chimera_solver")
chimeraSolver = chimera_solver_module.CreateSolver()



''' ##########################################  '''
''' SETUP for BACKGROUND  ''' 
''' ##########################################  '''
background_parameter_file = open(backGroundParameterFile,'r')
backGroundProjectParameters = Parameters( background_parameter_file.read())

## Fluid model part definition 
backGround_model_part = ModelPart(backGroundProjectParameters["problem_data"]["model_part_name"].GetString())
backGround_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, backGroundProjectParameters["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available
backGroundModel = {backGroundProjectParameters["problem_data"]["model_part_name"].GetString() : backGround_model_part}

## Solver construction
backGround_solver_module = __import__(backGroundProjectParameters["solver_settings"]["solver_type"].GetString())
back_ground_solver = backGround_solver_module.CreateSolver(backGround_model_part, backGroundProjectParameters["solver_settings"])
################ 
chimeraSolver.AddBackgroundModelPart(backGround_model_part, back_ground_solver)
#####################################
back_ground_solver.AddVariables()
## Read the model - note that SetBufferSize is done here
back_ground_solver.ImportModelPart()
## Add AddDofs
back_ground_solver.AddDofs()


## Initialize GiD  I/O
from gid_output_process import GiDOutputProcess
background_gid_output = GiDOutputProcess(back_ground_solver.GetComputeModelPart(),
                              backGroundProjectParameters["problem_data"]["problem_name"].GetString() ,
                              backGroundProjectParameters["output_configuration"])
background_gid_output.ExecuteInitialize()
##here all of the allocation of the strategies etc is done
back_ground_solver.Initialize()


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
print('################ :: setup of back ground model finished !')


''' ##########################################  '''
''' SETUP for PATCHES  ''' 
''' ##########################################  '''
for patchName, patchParameterFile in patchNameParameterFileMap.items():
    parameter_file = open(patchParameterFile,'r')
    ProjectParameters = Parameters( parameter_file.read())
    patchNameParameterMap[patchName] = ProjectParameters

for patchName, patchProjectParameters in patchNameParameterMap.items():
    patchModelPart = ModelPart(patchProjectParameters["problem_data"]["model_part_name"].GetString())
    patchModelPart.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())
    patchNameModelPartMap[patchName] = patchModelPart

for patchName, patchProjectParameters in patchNameParameterMap.items():
    patchModelPart = patchNameModelPartMap[patchName]
    patchModel = {patchProjectParameters["problem_data"]["model_part_name"].GetString() : patchModelPart}
    patchNameModelMap[patchName] = patchModel
    
for patchName, patchProjectParameters in patchNameParameterMap.items():
    patchModelPart = patchNameModelPartMap[patchName]    
    patch_solver_module = __import__(patchProjectParameters["solver_settings"]["solver_type"].GetString())
    patch_solver = patch_solver_module.CreateSolver(patchModelPart, patchProjectParameters["solver_settings"])
    chimeraSolver.AddChimeraPatchModelPart(patchName, patchModelPart, patch_solver)
    patch_solver.AddVariables()
    ## Read the model - note that SetBufferSize is done here
    patch_solver.ImportModelPart()
    ## Add AddDofs
    patch_solver.AddDofs()
    patchNameSolverMap[patchName] = patch_solver


for patchName, patchProjectParameters in patchNameParameterMap.items():
    patch_solver = patchNameSolverMap[patchName]    
    patch_gid_output = GiDOutputProcess(patch_solver.GetComputeModelPart(),
                              patchProjectParameters["problem_data"]["problem_name"].GetString() ,
                              patchProjectParameters["output_configuration"])
    patch_gid_output.ExecuteInitialize()
    patch_solver.Initialize()
    patchNameGidOutPutMap[patchName] = patch_gid_output
    
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


for patchName, patchProjectParameters in patchNameParameterMap.items():
    patchModel = patchNameModelMap[patchName]    
    patch_list_of_processes = process_factory.KratosProcessFactory(patchModel).ConstructListOfProcesses( patchProjectParameters["initial_conditions_process_list"] )
    patch_list_of_processes += process_factory.KratosProcessFactory(patchModel).ConstructListOfProcesses( patchProjectParameters["boundary_conditions_process_list"] )
    patch_list_of_processes += process_factory.KratosProcessFactory(patchModel).ConstructListOfProcesses( patchProjectParameters["gravity"] )
    ## Processes initialization
    for process in patch_list_of_processes:
        process.ExecuteInitialize()
    patchNameListProcessesMap[patchName] = patch_list_of_processes

print('################ :: setup of chimera patches finished !')                


''' ##########################################  '''
''' Helper Functions  ''' 
''' ##########################################  '''

def WriteInterfaces(time, testMPOuter, testMPInner):
    from gid_output import GiDOutput
    interfaceOut = GiDOutput("InterfaceResidualOuter"+str(time),True,"Ascii","Single",False,True);
    interfaceIn = GiDOutput("InterfaceResidualInner"+str(time),True,"Ascii","Single",False,True);
    interfaceOut.initialize_results(testMPOuter);
    interfaceIn.initialize_results(testMPInner)
    interfaceOut.write_results(time, testMPOuter, ["DISTANCE"], ["VELOCITY"])
    interfaceIn.write_results(time, testMPInner, ["DISTANCE"], ["VELOCITY"])


### The patch is moved in a sinusoidal way using this function.
### Period of oscillation is 0.10 second
### Amplitude of oscillation is 1.5
def MovePatchMesh(patch_model_part, time):
    perd = 0.10;
    amp = 0.2;
    fac = amp*math.sin(perd*np.pi*time)
    for node in patch_model_part.Nodes:
        node.X = node.X + fac


''' ##########################################  '''
''' CHIMERA THINGS  ''' 
''' ##########################################  '''
'''chimeraSolver.Initialize('Inlet3D_interface')
testMPOuter, testMPInner = chimeraSolver.UpdateOverlapBoundaries()

WriteInterfaces(0.0, testMPOuter, testMPInner)
testMPInner = None
testMPOuter = None'''


''' ##########################################  '''
''' SOLUTION PROCEDURE  ''' 
''' ##########################################  '''
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

print('#######################')
print('Starting Time Step Loop')
print('#######################')

while(time <= end_time):
    time = time + Dt
    step = step + 1
    print("STEP = ", step)
    print("TIME = ", time)
    backGround_model_part.CloneTimeStep(time)
    for patchName, patch_model_part in patchNameModelPartMap.items():
        patch_model_part.CloneTimeStep(time)
    
    
    #chimeraSolver.GetPressureOnChimeraParts()
    testMPOuter, testMPInner = chimeraSolver.UpdateOverlapBoundaries()
    WriteInterfaces(time, testMPOuter, testMPInner)
    testMPInner = None
    testMPOuter = None
    
    if(step >= 3):
        for process in background_list_of_processes:
            process.ExecuteInitializeSolutionStep()            
        for patchName, patch_list_of_processes in patchNameListProcessesMap.items():
            for process in patch_list_of_processes:
                process.ExecuteInitializeSolutionStep()
        
        background_gid_output.ExecuteInitializeSolutionStep()
        for patchName, patch_gid_output in patchNameGidOutPutMap.items():
            patch_gid_output.ExecuteInitializeSolutionStep()
                
            
        for patchName, patch_model_part in patchNameModelPartMap.items():
            MovePatchMesh(patch_model_part, time)
            
        ##### Chimera Solver
        chimeraSolver.Solve()    
        
        for process in background_list_of_processes:
            process.ExecuteFinalizeSolutionStep()            
        for patchName, patch_list_of_processes in patchNameListProcessesMap.items():
            for process in patch_list_of_processes:
                process.ExecuteFinalizeSolutionStep()        

        background_gid_output.ExecuteFinalizeSolutionStep()        
        for patchName, patch_gid_output in patchNameGidOutPutMap.items():
            patch_gid_output.ExecuteFinalizeSolutionStep()
            

        #TODO: decide if it shall be done only when output is processed or not
        for process in background_list_of_processes:
            process.ExecuteBeforeOutputStep()        
        for patchName, patch_list_of_processes in patchNameListProcessesMap.items():
            for process in patch_list_of_processes:
                process.ExecuteBeforeOutputStep()

    
        if background_gid_output.IsOutputStep():
            background_gid_output.PrintOutput()
        for patchName, patch_gid_output in patchNameGidOutPutMap.items():
            if patch_gid_output.IsOutputStep():
                patch_gid_output.PrintOutput()
        
        for process in background_list_of_processes:
            process.ExecuteAfterOutputStep()        
        for patchName, patch_list_of_processes in patchNameListProcessesMap.items():
            for process in patch_list_of_processes:
                process.ExecuteAfterOutputStep()        

        out = out + Dt


for process in background_list_of_processes:
    process.ExecuteFinalize()        
for patchName, patch_list_of_processes in patchNameListProcessesMap.items():
    for process in patch_list_of_processes:
        process.ExecuteFinalize()       


background_gid_output.ExecuteFinalize()        
for patchName, patch_gid_output in patchNameGidOutPutMap.items():
    patch_gid_output.ExecuteFinalize()    
