
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
sys.path.append(kratos_benchmarking_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_StructuralApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)
from KratosStructuralApplication import *


import benchmarking

def FindNode(node_list,x,y,z):
    for node in node_list:
        if ((node.X - x) ** 2 + (node.Y - y) ** 2 + (node.Z - z) ** 2 < 0.0000001):
            print  node
            return node
    
def BenchmarkCheck(time, node1):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_X), "Node 1 Desplacement_x", 1.0)
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_Y), "Node 1 Desplacement_y", 1.0)


## from now on the order is not anymore crucial
#defining a model part
model_part = ModelPart("FluidPart");  
model_part.AddNodalSolutionStepVariable(PK2_STRESS_TENSOR);
model_part.AddNodalSolutionStepVariable(BODY_FORCE);

#adding of Variables to Model Part should be here 
import structural_solver_dynamic
structural_solver_dynamic.AddVariables(model_part)



#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("Pendulo_Kratos",gid_mode,multifile,deformed_mesh_flag, write_conditions)
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
structural_solver_dynamic.AddDofs(model_part)


model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
model_part.Properties[2].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
print "Linear elastic model selected"

#creating a fluid solver object
solver = structural_solver_dynamic.DynamicStructuralSolver(model_part,domain_size)





node_1 = FindNode(model_part.Nodes, 1.53624,   1.53022,   0.00000)

solver.Initialize()

Dt = 0.1
nsteps = 100

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())


for step in range(0,nsteps):
    print "line49"

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 0):
        solver.Solve()
        if (benchmarking.InBuildReferenceMode()):
	 BenchmarkCheck(time, node_1)
        else:
         BenchmarkCheck(time, node_1)

        #print the results
        gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
        gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(ACCELERATION,model_part.Nodes,time,0)
#        gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time,domain_size)
    

gid_io.FinalizeResults()
print "COMPLETED ANALYSIS"
    
        






