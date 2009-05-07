##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3


##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path


kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
kratos_python_scripts_path = '../../../../applications/structural_application/python_scripts'
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking

import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_python_scripts_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_StructuralApplication = True
applications_interface.Import_KratosExternalSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)
from KratosStructuralApplication import *
print kernel

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosExternalSolversApplication import *

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here 
import structural_solver_static

structural_solver_static.AddVariables(model_part)

###################################################################
###################################################################

model_part.AddNodalSolutionStepVariable(FORCE);
model_part.AddNodalSolutionStepVariable(DENSITY);
model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
model_part.AddNodalSolutionStepVariable(REACTION);

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
#gid_mode = GiDPostMode.GiD_PostAscii
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

gid_io = GidIO("arc_length_des",gid_mode,multifile,deformed_mesh_flag, write_conditions)
gid_io.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()

print model_part

#writing the mesh
#gid_io.WriteMesh(model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
structural_solver_static.AddDofs(model_part)

#creating a fluid solver object
solver = structural_solver_static.StaticStructuralSolver(model_part,domain_size)

solver.structure_linear_solver =  SuperLUSolver()

solver.Initialize()
(solver).SetEchoLevel(2);

print "Forcing Node List"
forcing_nodes_list = [];

for node in model_part.Nodes:
    if(node.Y > 0.249 and node.Y < 0.2501):
        forcing_nodes_list.append(node)

for node in forcing_nodes_list:
          print node.Id , " ", node.X, " ",node.Y

model_part.Nodes[1].SetSolutionStepValue(FORCE_Y,0,0.00);


def IncreaseDisplacement(forcing_nodes_list,disp):
    for node in forcing_nodes_list:
        node.SetSolutionStepValue(DISPLACEMENT_Y,0,disp)
        
 

Dt         = 1;
nsteps     = 60;
Inc        = 0.01;
disp       = 0.00;

print("initializing results")

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())


for inner_step in range(1,nsteps):

    time = Dt*inner_step
    model_part.CloneTimeStep(time)

    print time
    disp = -Inc*inner_step;
    IncreaseDisplacement(forcing_nodes_list,disp);
    print disp; 
    solver.Solve() 

                   
    #print the results
    print " Writing Result in Gid"
    gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(FORCE,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
    gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time)
    gid_io.PrintOnGaussPoints(DAMAGE,model_part,time)

gid_io.FinalizeResults()
print "COMPLETED ANALYSIS"   
