##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2
import math

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
applications_interface.Import_StructuralApplication = True
applications_interface.Import_KratosExternalSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)


import benchmarking

## from now on the order is not anymore crucial
##################################################################
##################################################################


def FindNode(node_list,x,y,z):
    for node in node_list:
        if ((node.X - x) ** 2 + (node.Y - y) ** 2 + (node.Z - z) ** 2 < 0.0000001):
            print  node
            return node
    
def BenchmarkCheck(time, node1, node2):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_X), "Node 1 Displacement_x", 1.0)
    benchmarking.Output(node2.GetSolutionStepValue(DISPLACEMENT_Y), "Node 2 Displacement_y", 1.0)

def AnalyticalResults(time, node1, node2):
    benchmarking.Output(time, "Time")
    benchmarking.Output(0.00006, "Node 1 Displacement_x", 1.0)
    benchmarking.Output(0.00024, "Node 2 Displacement_y", 1.0)




from KratosStructuralApplication import *
print kernel

from KratosExternalSolversApplication import *

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here 
#import structural_solver_dynamic
#import structural_solver_dynamic_superlu
#import structural_solver_relaxation
import structural_solver_static

structural_solver_static.AddVariables(model_part)
#structural_solver_relaxation.AddVariables(model_part)
#structural_solver_dynamic.AddVariables(model_part)
#structural_solver_dynamic_superlu.AddVariables(model_part)

###################################################################
###################################################################

model_part.AddNodalSolutionStepVariable(FORCE);
model_part.AddNodalSolutionStepVariable(DENSITY);
model_part.AddNodalSolutionStepVariable(DISPLACEMENT);


#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
#gid_mode = GiDPostMode.GiD_PostAscii
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

gid_io = GidIO("Patch_Test_Total_Lagrangian_4N",gid_mode,multifile,deformed_mesh_flag, write_conditions)
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
#structural_solver_dynamic.AddDofs(model_part)
#structural_solver_dynamic_superlu.AddDofs(model_part)
#structural_solver_relaxation.AddDofs(model_part)
structural_solver_static.AddDofs(model_part)

#creating a fluid solver object
#solver = structural_solver_relaxation.RelaxationStructuralSolver(model_part,domain_size)
#solver = structural_solver_dynamic.DynamicStructuralSolver(model_part,domain_size)
#solver = structural_solver_dynamic_superlu.DynamicStructuralSolver(model_part,domain_size)
solver = structural_solver_static.StaticStructuralSolver(model_part,domain_size)

solver.structure_linear_solver =  SuperLUSolver()


model_part.Properties[1].SetValue(CONSTITUTIVE_LAW,  Isotropic2D())
print "Isotropic model selected"


node_1 = FindNode(model_part.Nodes, 0.00, 0.12, 0.00)
node_2 = FindNode(model_part.Nodes, 0.24, 0.12, 0.00)


solver.Initialize()
(solver).SetEchoLevel(2);


Dt         = 1;
nsteps     = 4;

print("initializing results")

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())


for step in range(0,nsteps):

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    if(step>2):
        solver.Solve()
	if (benchmarking.InBuildReferenceMode()):
	  AnalyticalResults(time, node_1, node_2)
	else:
	  BenchmarkCheck(time, node_1, node_2)




 
	#print the results
        print " Writing Result in Gid"
        gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(FORCE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
        gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time)


gid_io.FinalizeResults()
print "COMPLETED ANALYSIS"   
