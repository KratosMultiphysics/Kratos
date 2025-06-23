from math import exp,sqrt
##################################################################
##################################################################
#setting the domain size for the problem to be solved


##################################################################
##################################################################

#including kratos path
from KratosMultiphysics import *    #we import the KRATOS  
from KratosMultiphysics.ConvectionDiffusionApplication import *        #and now our application. note that we can import as many as we need to solve our specific problem 

#import eulerian_convection_diffusion_solver as convection_diffusion_solver
#import pfem2_convection_diffusion_solver as convection_diffusion_solver
import bfecc_convection_diffusion_solver as convection_diffusion_solver

import ProjectParameters

ConvDiffSettings = ProjectParameters.ConvectionSolverSettings

domain_size = ProjectParameters.domain_size 

#defining a model part
model_part = ModelPart("DummyPart");  

#THERE ARE TWO WAYS TO DEFINE THE CONVECTION DIFFUSION SETTINGS.
#OPTION 1 is using the ProjectParameters , where the parameters are defined like unknown_variable="TEMPERATURE" and so on:
#this option is meant for GUIs, although not so transparent to the user
ConvDiffSettings = ProjectParameters.ConvectionSolverSettings
convection_diffusion_solver.AddVariables(model_part, ConvDiffSettings)

#OPTION 2 is aimed at developers. The CONVECTION_DIFFUSION_SETTINGS of the model_part are defined by the user in this file
#then the AddVariables will read the values directly from the model part: 
#thermal_settings = ConvectionDiffusionSettings()
#thermal_settings.SetDiffusionVariable(CONDUCTIVITY)
#thermal_settings.SetUnknownVariable(TEMPERATURE)
#thermal_settings.SetVelocityVariable(VELOCITY)
#(model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, thermal_settings)
#eulerian_convection_diffusion_solver.AddVariables(model_part) # no settings here! (read from the model_part)


#now we proceed to use the GID interface (both to import the information inside the .mdpa file and later print the results in a file  
gid_mode = GiDPostMode.GiD_PostBinary  #we import the python file that includes the commands that we need  
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("PecletTest",gid_mode,multifile,deformed_mesh_flag,write_conditions) #output will be the name of the output files

model_part_io = ModelPartIO("PecletTest")             # we set the name of the .mdpa file  
model_part_io.ReadModelPart(model_part)                                                 # we load the info from the .mdpa 

conductivity=0.1
velocity=1.0

#the ConvDiffSettings can be passed, but are ignored in the AddDofs. (only added to match 
#After adding the variables, only the model_part.convdiffsettings is used.
convection_diffusion_solver.AddDofs(model_part)

#BOUNDARY CONDITIONS MUST ALSO BE SET BY THE USER!
for node in model_part.Nodes:
    node.SetSolutionStepValue(VELOCITY_X,velocity)
    node.SetSolutionStepValue(CONDUCTIVITY,conductivity)
    node.SetSolutionStepValue(TEMPERATURE,node.X)
    if node.X<0.0001:
        node.Fix(TEMPERATURE);
    if node.X>0.9999:
        node.Fix(TEMPERATURE);


# we create a mesh (to be used in the postprocess)
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
print(model_part)

#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)  
model_part.SetBufferSize(2)

convection_solver = convection_diffusion_solver.CreateSolver(model_part, ProjectParameters)
convection_solver.Initialize()

gid_io.InitializeResults(mesh_name,(model_part).GetMesh()) 
print ("hola2")
time =0.0
delta_time=0.1
final_time=10.0
step = 0
while time<final_time:
   print ("current step",step)
   step += 1
   time += delta_time
   model_part.CloneTimeStep(time)
   if step>3:
        convection_solver.Solve()
   gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,time,0)
   

gid_io.FinalizeResults()

error = 0.0
nnodes = 0.0
for node in model_part.Nodes:
    Temp = node.GetSolutionStepValue(TEMPERATURE)
    aTemp = (1.0-exp(velocity*node.X/conductivity))/(1.0-exp(velocity/conductivity))
    error += (Temp-aTemp)**2
    print("the error at x=",node.X," is ", Temp-aTemp )
    nnodes +=1.0

error/= nnodes
error = sqrt(error)
print("the error is ",error)

