##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

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
sys.path.append(kratos_benchmarking_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_ExternalSolversApplication = True
applications_interface.Import_StructuralApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

import benchmarking

def FindNode(node_list,x,y,z):
    for node in node_list:
        if ((node.X - x) ** 2 + (node.Y - y) ** 2 + (node.Z - z) ** 2 < 0.0000001):
            print  node
            return node
    
def BenchmarkCheck(time, node1, node2, node3, node4):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_X), "Node 1 Displacement_x", 0.00001)
    benchmarking.Output(node2.GetSolutionStepValue(DISPLACEMENT_Y), "Node 2 Displacement_y", 0.00001)
    benchmarking.Output(node3.GetSolutionStepValue(REACTION_X), "Node 3 Reaction_x", 0.00001)
    benchmarking.Output(node4.GetSolutionStepValue(REACTION_Y), "Node 4 Reaction_y", 0.00001)



#def AnalyticalResults(time, node1, node2,node3, node4):
    #benchmarking.Output(time, "Time")
    #benchmarking.Output(-0.221921365586,  "Node 1 Displacement_x", 0.00001)
    #benchmarking.Output(-0.0361068223759, "Node 2 Displacement_y", 0.00001)
    #benchmarking.Output( 51.6844785228,   "Node 3 Reaction_x", 0.00001)
    #benchmarking.Output( -123.134969306,   "Node 4 Reaction_y", 0.00001)





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
model_part.AddNodalSolutionStepVariable(BODY_FORCE);

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
 
gid_io = GidIO("cantilever2d",gid_mode,multifile,deformed_mesh_flag, write_conditions)
gid_io.ReadModelPart(model_part)
print "writing mesh..."
print model_part
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
print "pippo"
print model_part
#print model_part.Properties

#writing the mesh
#gid_io.WriteMesh(model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
structural_solver_static.AddDofs(model_part)

#creating a fluid solver object
solver = structural_solver_static.StaticStructuralSolver(model_part,domain_size)
##pILUPrecond = ILU0Preconditioner() 
##solver.structure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
solver.structure_linear_solver = SuperLUSolver()
model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
print "Linear elastic model selected"

node_1 = FindNode(model_part.Nodes, 0.05, 1.00, 0.00)
node_2 = FindNode(model_part.Nodes, 0.00, 1.00, 0.00)
node_3 = FindNode(model_part.Nodes, 0.00, 0.00, 0.00)
node_4 = FindNode(model_part.Nodes, 0.05, 0.00, 0.00)

solver.Initialize()
(solver).SetEchoLevel(2);

gravity     = -9.80665
for node in model_part.Nodes:
    node.SetSolutionStepValue(BODY_FORCE_X,0,0.0);
    node.SetSolutionStepValue(BODY_FORCE_Y,0,gravity);
    node.SetSolutionStepValue(BODY_FORCE_Z,0,0.0); 

Dt = 0.001
nsteps =9
print("initializing results")
gid_io.InitializeResults(mesh_name,(model_part).GetMesh())
for step in range(0,nsteps):

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
       solver.Solve()
       if (benchmarking.InBuildReferenceMode()):
         #AnalyticalResults(time, node_1, node_2, node_3, node_4)
	 BenchmarkCheck(time, node_1, node_2, node_3, node_4)
       else:
         BenchmarkCheck(time, node_1, node_2, node_3, node_4)

    #print the results
    print "a"
    gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
gid_io.FinalizeResults()
print "finito"


          
        

