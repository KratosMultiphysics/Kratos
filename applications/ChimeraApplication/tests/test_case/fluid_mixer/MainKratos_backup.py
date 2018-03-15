from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
from KratosMultiphysics import *

from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

from KratosMultiphysics.ChimeraApplication import *
import csv 
import math
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
patch1 = main_model_part.GetSubModelPart("GENERIC_patch1")
patch2 = main_model_part.GetSubModelPart("GENERIC_patch2")
patchBoundary1 = main_model_part.GetSubModelPart("GENERIC_patchBoundary1")
patchBoundary2 = main_model_part.GetSubModelPart("GENERIC_patchBoundary2")
structure1 = main_model_part.GetSubModelPart("NoSlip2D_structure1")
structure2 = main_model_part.GetSubModelPart("NoSlip2D_structure2")
#background.AddNodalSolutionStepVariable(DISTANCE)
#patch.AddNodalSolutionStepVariable(DISTANCE)

## Initialize GiD  I/O
from gid_output_process import GiDOutputProcess
gid_output_backGround = GiDOutputProcess(background,
                              "BackgroundSolution" ,
                              ProjectParameters["output_configuration"])

gid_output_patch1 = GiDOutputProcess(patch1,
                              "PatchSolution" ,
                              ProjectParameters["output_configuration"])
gid_output_patch2 = GiDOutputProcess(patch2,
                              "PatchSolution" ,
                              ProjectParameters["output_configuration"])

gid_output_backGround.ExecuteInitialize()
gid_output_patch1.ExecuteInitialize()
gid_output_patch2.ExecuteInitialize()

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
vtkOutput_patch1 = VtkOutput(patch1,"nnn",ProjectParameters["output_configuration"])
vtkOutput_patch2 = VtkOutput(patch2,"nnn",ProjectParameters["output_configuration"])
vtkOutput_background = VtkOutput(background,"nnn",ProjectParameters["output_configuration"])

#TODO: think if there is a better way to do this
#fluid_model_part = solver.GetComputeModelPart()

## Stepping and time settings
Dt = ProjectParameters["problem_data"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

time = 0.0
step = 0
out = 0.0

gid_output_backGround.ExecuteBeforeSolutionLoop()
gid_output_patch1.ExecuteBeforeSolutionLoop()
gid_output_patch2.ExecuteBeforeSolutionLoop()

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
    
    
file_name_patch = 'patch.csv' 
file_name_background = 'background.csv'  


def PrintDrag(model_part, time,Dt):
  
  filename = 'drag.csv'
  
  if time <= 3*Dt:
    with open (filename,'w') as file:
      file.write('Time,RX,RY,RZ,MX,MY,MZ\n')
  
   
  rx = 0.0
  ry = 0.0
  rz = 0.0

  mx = 0.0
  my = 0.0
  mz = 0.0
  
  centre = [0,0,0]
  n_nodes = len(model_part.Nodes)
  
  for node in model_part.Nodes:    
    centre[0] += node.X
    centre[1] += node.Y
    centre[2] += node.Z
    
  centre[0] = centre[0]/n_nodes
  centre[1] = centre[1]/n_nodes
  centre[2] = centre[2]/n_nodes
  
  


  for node in model_part.Nodes:
      reaction = node.GetSolutionStepValue(REACTION, 0)
      rx += reaction[0]
      ry += reaction[1]
      rz += reaction[2]
      

      x = node.X - centre[0]
      y = node.Y - centre[1]
      z = node.Z - centre[2]
      mx -= y * reaction[2] - z * reaction[1]
      my -= z * reaction[0] - x * reaction[2]
      mz -= x * reaction[1] - y * reaction[0]

  with open (filename,'a') as file:
    file.write(str(time) + "," + str(rx) + "," + str(ry) + "," + str(rz) + "," + str(mx) + "," + str(my) + "," + str(mz) + "\n")
    
def RevolvePatch(model_part,Dt,omega,Rc,Ro):
  

  theta =  omega*Dt


  costheta = math.cos(theta)
  sintheta = math.sin(theta)

  rox = Ro[0]-Rc[0]
  roy = Ro[1]-Rc[1]
  
  dx = (costheta-1)*rox - sintheta*roy
  dy = sintheta*rox + (costheta-1)*roy
  
  print('dx',dx)
  print('dy',dy)

  Ro[0] = Ro[0]+dx
  Ro[1] = Ro[1]+dy

  x = Ro[0]-Rc[0]
  y = Ro[1]-Rc[1]



  vel_x = -omega*y
  vel_y = omega*x
  
  for node in model_part.Nodes:
    node.X = node.X + dx
    node.Y = node.Y + dy
    node.X0 = node.X0 + dx
    node.Y0 = node.Y0 + dy
    node.SetSolutionStepValue(MESH_VELOCITY_X,0,vel_x)
    node.SetSolutionStepValue(MESH_VELOCITY_Y,0,vel_y)
    
   
  
  
    

ChimeraParamaters1 = Parameters(R"""{
                                "process_name":"apply_chimera_process",
                                "background": {
                                                "model_part_name":"GENERIC_background",
                                                 "pressure_coupling":"all",
                                                  "type" : "nearest_element",
                                                  "IsWeak" : true
                                                },
                                "patch" :       {
                                                  "model_part_name":"GENERIC_patch1",
                                                  "pressure_coupling" : "all",
                                                  "type" : "nearest_element",
                                                  "IsWeak" : true
                                                },
                                "pressure_coupling_node" : 0.0,
                                "patch_boundary_model_part_name":"GENERIC_patchBoundary1",
                                "overlap_distance":0.05
                                }""")

ChimeraParamaters2 = Parameters(R"""{
                                "process_name":"apply_chimera_process",
                                "background": {
                                                "model_part_name":"GENERIC_background",
                                                 "pressure_coupling":"all",
                                                  "type" : "nearest_element",
                                                  "IsWeak" : true
                                                },
                                "patch" :       {
                                                  "model_part_name":"GENERIC_patch2",
                                                  "pressure_coupling" : "all",
                                                  "type" : "nearest_element",
                                                  "IsWeak" : true
                                                },
                                "pressure_coupling_node" : 0.0,
                                "patch_boundary_model_part_name":"GENERIC_patchBoundary2",
                                "overlap_distance":0.05
                                
                                }""")

ChimeraProcess1 = ApplyChimeraProcess2d(main_model_part,ChimeraParamaters1)
ChimeraProcess2 = ApplyChimeraProcess2d(main_model_part,ChimeraParamaters2)


RotateRegionParameters1 = Parameters(R"""{
                "movement_name":"default",
                "sub_model_part_name":"GENERIC_patch1",
                "sub_model_part_boundary_name":"GENERIC_patchBoundary1",          
                "center_of_rotation":[-1,0,0],
                "angular_velocity_radians":6.283184,
                "axis_of_rotation":[0,0,1],
                "is_ale" : true
            }""")
RotateRegionParameters2 = Parameters(R"""{
                "movement_name":"default",
                "sub_model_part_name":"GENERIC_patch2",
                "sub_model_part_boundary_name":"GENERIC_patchBoundary2",          
                "center_of_rotation":[1,0,0],
                "angular_velocity_radians":-6.283184,
                "axis_of_rotation":[0,0,1],
                "is_ale" : true
            }""")
RotateRegionProcess1 = RotateRegionProcess(main_model_part,RotateRegionParameters1)
RotateRegionProcess2 = RotateRegionProcess(main_model_part,RotateRegionParameters2)
patch_centre1 = [-1,0.0,0.0]
patch_centre2 = [1,0.0,0.0]
orgin = [0.0,0.0,0.0]

while(time <= end_time):

    time = time + Dt
    step = step + 1 
    main_model_part.CloneTimeStep(time)


    print("STEP = ", step)
    print("TIME = ", time)
    
    main_model_part.ProcessInfo[STEP]=step

    
    if(step >= 3):
        #print("Nr of processes",len(list_of_processes))
        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()
            #print("process ", process)
        main_model_part.GetNode(1).Fix(PRESSURE)
        main_model_part.GetNode(1).SetSolutionStepValue(PRESSURE,0,0)
        RevolvePatch(patch1,Dt,3.141592,orgin,patch_centre1)
        RevolvePatch(patch2,Dt,3.141592,orgin,patch_centre2)
        

        RotateRegionProcess1.SetCentreOfRotation(patch_centre1[0],patch_centre1[1],patch_centre1[2])
        RotateRegionProcess2.SetCentreOfRotation(patch_centre2[0],patch_centre2[1],patch_centre2[2])

        RotateRegionProcess1.ExecuteInitializeSolutionStep()  
        RotateRegionProcess2.ExecuteInitializeSolutionStep() 
        ChimeraProcess1.ExecuteInitializeSolutionStep()
        ChimeraProcess2.ExecuteInitializeSolutionStep()

            
        gid_output_backGround.ExecuteInitializeSolutionStep()
        gid_output_patch1.ExecuteInitializeSolutionStep()
        gid_output_patch2.ExecuteInitializeSolutionStep()
        
        
        
        
        solver.Solve()
        PrintDrag(structure1, time,Dt)
        PrintDrag(structure2, time,Dt)

        #print("exiting the process loop")
        
        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()
            
        gid_output_backGround.ExecuteFinalizeSolutionStep()
        gid_output_patch1.ExecuteFinalizeSolutionStep()
        gid_output_patch2.ExecuteFinalizeSolutionStep()

        #TODO: decide if it shall be done only when output is processed or not
        for process in list_of_processes:
          process.ExecuteBeforeOutputStep()
    
        if gid_output_backGround.IsOutputStep():
          gid_output_backGround.PrintOutput()
          vtkOutput_background.PrintOutput()
            
        if gid_output_patch1.IsOutputStep():
          gid_output_patch1.PrintOutput()
          vtkOutput_patch1.PrintOutput()  
            
        if gid_output_patch2.IsOutputStep():
          gid_output_patch2.PrintOutput()
          vtkOutput_patch2.PrintOutput() 
        for process in list_of_processes:
          process.ExecuteAfterOutputStep()
            

      
        RotateRegionProcess1.ExecuteFinalizeSolutionStep()
        RotateRegionProcess2.ExecuteFinalizeSolutionStep()
        ChimeraProcess1.ExecuteFinalizeSolutionStep()
        ChimeraProcess2.ExecuteFinalizeSolutionStep()

        out = out + Dt

for process in list_of_processes:
    process.ExecuteFinalize()
    
gid_output_backGround.ExecuteFinalize()
gid_output_patch1.ExecuteFinalize()
gid_output_patch2.ExecuteFinalize()














