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
applications_interface.Import_KratosMKLSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)
from KratosStructuralApplication import *
print kernel
## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosMKLSolversApplication import *

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here
import structural_solver_static
structural_solver_static.AddVariables(model_part)

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
 
gid_io = GidIO("cantilever2d",gid_mode,multifile,deformed_mesh_flag, write_conditions)
gid_io.ReadModelPart(model_part)
print "writing mesh..."
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
print model_part
#print model_part.Properties

#writing the mesh
#gid_io.WriteMesh(model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
structural_solver_static.AddDofs(model_part)

model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
print "Linear elastic model selected"

#creating a fluid solver object
solver = structural_solver_static.StaticStructuralSolver(model_part,domain_size)
#pILUPrecond = ILU0Preconditioner() 
#solver.structure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)



solver.structure_linear_solver =  MKLPardisoSolver()


solver.Initialize()
(solver).SetEchoLevel(2);

Dt = 0.001
nsteps = 5
print("initializing results")
gid_io.InitializeResults(mesh_name,(model_part).GetMesh())
for step in range(0,nsteps):

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    solver.Solve()

    #print the results
    print "a"
    gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
    gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time)
gid_io.FinalizeResults()
print "finito"


          
        

