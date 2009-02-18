##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs/' ##kratos_root/libs
kratos_applications_path = '../../../../applications/' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_MeshingApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

#loading meshing application
sys.path.append(kratos_applications_path + 'meshing_application/python_scripts') 
from KratosMeshingApplication import *
meshing_application = KratosMeshingApplication()
kernel.AddApplication(meshing_application)

applications_interface.ImportApplications(kernel, kratos_applications_path)

#loading meshing application
kernel.InitializeApplication(meshing_application);

## from now on the order is not anymore crucial
##################################################################
##################################################################

import math

#defining a model part
model_part = ModelPart("FluidPart");  

#adding variables to be used
model_part.AddNodalSolutionStepVariable(NODAL_H)
model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
model_part.AddNodalSolutionStepVariable(TEMPERATURE)
model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
model_part.AddNodalSolutionStepVariable(IS_FLUID)


#reading the model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("adaptive_mesher2d",gid_mode_flag,use_multifile,deformed_print_flag,write_conditions)
gid_io.ReadModelPart(model_part)

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#giving a initial "temperature" distribtion
for node in model_part.Nodes:
    temp = node.X**2 + node.Y**2 + node.Z**2
    node.SetSolutionStepValue(TEMPERATURE,0,temp);
    node.SetSolutionStepValue(DISPLACEMENT_X,0,node.X);
    
#defining the mesher
Mesher = TriGenPFEMModeler()
#Mesher = TriGenModeler()

for node in model_part.Nodes:
    node.SetSolutionStepValue(NODAL_H,0, 0.1)
    node.SetSolutionStepValue(IS_FLUID,0, 1)
    
node_erase_process = NodeEraseProcess(model_part);

neigh_finder = FindNodalNeighboursProcess(model_part,9,18)
neigh_finder.Execute()
    
Dt = 0.01
nsteps = 12
time = 0.0
i=0
for step in range(0,nsteps):

    time = time + Dt

    print "time = ", time, " new_Dt= ",Dt," step = ", step
    
    model_part.CloneTimeStep(time)
    
    print "TimeStep cloned!!!"

    if(step >= 3) :
        #erasing all of the elements
        (model_part.Elements).clear();
        (model_part.Conditions).clear();
        
        print "Arrays are cleared!!!"
        
        #remeshing
        alpha_shape = 1.8

        
        for node in model_part.Nodes:
            if ( ((node.X-1.0-i)*(node.X-1.0-i)) + ((node.Y-0.5)*(node.Y-0.5)) < 0.06):
                node.SetSolutionStepValue(NODAL_H,0, 0.01)
            else:
                node.SetSolutionStepValue(NODAL_H,0, 0.1)
        
        i=i+1;    
        
        Mesher.ReGenerateMesh("TestElement2D","Condition2D", model_part, node_erase_process, alpha_shape)
    
        print "meshing is performed"
        
        #recalculating neighbours 
        (neigh_finder).Execute();
        
        #printing output
        file_name = 'time_' + str(time)
##        gid_io.ChangeOutputName(file_name,GiDPostMode.GiD_PostBinary);
        mesh_name = time

        gid_io.InitializeMesh( time );
        gid_io.WriteNodeMesh((model_part).GetMesh());
        gid_io.WriteMesh((model_part).GetMesh());
        gid_io.FinalizeMesh();

        gid_io.InitializeResults(time, (model_part).GetMesh());
        gid_io.WriteNodalResults(TEMPERATURE, model_part.Nodes, time, 0);
        gid_io.WriteNodalResults(NODAL_H, model_part.Nodes, time, 0);
        gid_io.WriteNodalResults(IS_BOUNDARY, model_part.Nodes, time, 0);
        gid_io.WriteNodalResults(DISPLACEMENT, model_part.Nodes, time, 0);

        gid_io.Flush()
        #gid_io.CloseResultFile();
        gid_io.FinalizeResults()

        print "step finished"

print "finito"



          
        

