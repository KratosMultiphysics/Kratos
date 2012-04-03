##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_path = '../../../../' ##kratos_root/
import sys
sys.path.append(kratos_path)

#importing Kratos main library
from KratosMultiphysics import *

#loading meshing application
sys.path.append(kratos_applications_path + 'meshing_application/python_scripts') 
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *



## from now on the order is not anymore crucial
##################################################################
##################################################################

# import math

# #defining a model part
# model_part = ModelPart("FluidPart");  

# #adding variables to be used
# model_part.AddNodalSolutionStepVariable(NODAL_H)
# model_part.AddNodalSolutionStepVariable(TEMPERATURE)
# model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
# model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
# model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
# model_part.AddNodalSolutionStepVariable(IS_FLUID)

# #reading the model
# gid_io = GidIO("dam3d_remesh",GiDPostMode.GiD_PostBinary)
# gid_io.ReadModelPart(model_part)
# print model_part

# #the buffer size should be set up here after the mesh is read for the first time
# model_part.SetBufferSize(3)

# #giving a initial "temperature" distribtion
# for node in model_part.Nodes:
#     temp = node.X**2 + node.Y**2 + node.Z**2
#     node.SetSolutionStepValue(TEMPERATURE,0,temp);

# #assigning a "nodal_h"
# for node in model_part.Nodes:
#     node.SetSolutionStepValue(NODAL_H,0,0.05);

# #defining the mesher
# Mesher = TetGenModeler()


# neigh_finder = FindNodalNeighboursProcess(model_part,20,30)
# neigh_finder.Execute()
    
# Dt = 0.005
# nsteps = 6
# time = 0.0

# for step in range(0,nsteps):

#     time = time + Dt

#     print "time = ", time, " new_Dt= ",Dt," step = ", step
    
#     model_part.CloneTimeStep(time)
    
#     print "TimeStep cloned!!!"
    
#     #erasing all of the elements
#     (model_part.Elements).clear();
#     (model_part.Conditions).clear();
    
#     print "Arrays are cleared!!!"
    
#     #remeshing
#     alpha_shape = 1.3

# ########## In this moment don't use this methdo !!!!!!!!!!!!!!
# #    Mesher.ReGenerateMeshElements(model_part,alpha_shape)

#     Mesher.ReGenerateMesh(model_part,alpha_shape)
    
#     print "mesing is performed"
    
#     #recalculating neighbours 
#     (neigh_finder).Execute();
    
#     #printing output
#     file_name = 'time_' + str(time)
#     gid_io.ChangeOutputName(file_name,GiDPostMode.GiD_PostBinary);
#     gid_io.WriteMesh((model_part).GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);
#     gid_io.WriteNodalResults(TEMPERATURE, model_part.Nodes, time, 0);
#     gid_io.Flush()



          
        

