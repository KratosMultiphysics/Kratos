from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

def Wait():
	input("Press Something")

#including kratos path
import sys
from KratosMultiphysics import *
#import KratosMultiphysics    #we import the KRATOS  
import KratosMultiphysics.MyLaplacianApplication as Poisson       #and now our application. note that we can import as many as we need to solve our specific problem 
import KratosMultiphysics.ExternalSolversApplication

#setting the domain size for the problem to be solved
domain_size = 2  # 2D problem  

#defining a model part
model_part = ModelPart("example");  #we create a model part  

import pure_diffusion_solver           #we import the python file that includes the commands that we need 
parameter_file = open("ProjectParameters.json",'r')  # we read the ProjectParameters.json file
custom_settings = Parameters(parameter_file.read())

pure_diffusion_solver = pure_diffusion_solver.CreateSolver(model_part, custom_settings)
pure_diffusion_solver.AddVariables()  #from the static_poisson_solver.py we call the function Addvariables so that the model part we have just created has the needed variables 

 # (note that our model part does not have nodes or elements yet) 

 #now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file  
gid_mode = GiDPostMode.GiD_PostAscii  #we import the python file that includes the commands that we need  
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("art4", gid_mode,multifile, deformed_mesh_flag, write_conditions)

model_part_io = ModelPartIO("example")          # we set the name of the .mdpa file  
model_part_io.ReadModelPart(model_part)         # we load the info from the .mdpa 

# we create a mesh for the postprocess  
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name )
gid_io.WriteMesh((model_part).GetMesh())
gid_io.FinalizeMesh()

#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)  
model_part.SetBufferSize(1)

# we add the DoFs  
pure_diffusion_solver.AddDofs()

# Initialize the solver
pure_diffusion_solver.Initialize()

print ("about to solve!")    
pure_diffusion_solver.Solve()
print ("Solved!")  

#and we print the results  
gid_io.InitializeResults(mesh_name,(model_part).GetMesh()) 
gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,0,0)
gid_io.FinalizeResults()

# since we have already calculated the temp, we can get the mean value
# first the constructor (it could have been called before)
calc_mean = Poisson.CalculateMeanTemperature(model_part)
# and we calculate!
calc_mean.Execute()