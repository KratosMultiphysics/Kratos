##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../../libs/' ##kratos_root/libs
kratos_applications_path = '../../../../../applications/' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_StructuralApplication = True
#applications_interface.Import_ExternalSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)
from KratosStructuralApplication import *

## from now on the order is not anymore crucial
##################################################################
##################################################################

#defining a model part 
model_part = ModelPart("FluidPart");  
 
#adding of Variables to Model Part should be here when the "very fix container will be ready" 
import structural_solver_dynamic
#import structural_solver_dynamic_superlu
structural_solver_dynamic.AddVariables(model_part) 
#structural_solver_dynamic_superlu.AddVariables(model_part)



#reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.SingleFile
deformed_print_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("cantilever3d",gid_mode_flag,use_multifile,deformed_print_flag,write_conditions)
gid_io.ReadModelPart(model_part)
print model_part 
print model_part.Properties 
 
#writing the mesh 
#gid_io.WriteMesh(model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary); 
 
#the buffer size should be set up here after the mesh is read for the first time 
model_part.SetBufferSize(2) 
 
#Adding DOFs
structural_solver_dynamic.AddDofs(model_part)
#structural_solver_dynamic_superlu.AddDofs(model_part)


model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )

#creating a fluid solver object
solver = structural_solver_dynamic.DynamicStructuralSolver(model_part,domain_size)
#testing superLU
#solver = structural_solver_dynamic_superlu.DynamicStructuralSolver(model_part, domain_size)

print "Linear elastic model selected"
print "Linear elastic model selected"
print "Linear elastic model selected"
print "Linear elastic model selected"
print "Linear elastic model selected"
print "Linear elastic model selected"
print "Linear elastic model selected"

print "Linear elastic model selected"
print "Linear elastic model selected"
print "Linear elastic model selected"
print "Linear elastic model selected"
print "Linear elastic model selected"
print "Linear elastic model selected"
print "Linear elastic model selected"
#gid_io.FinalizeMesh()

print "OLALALA"
solver.Initialize()

Dt = 0.02
nsteps = 8

print "OLALALA"

gid_io.InitializeMesh( 0.0 )
gid_io.WriteMesh( model_part.GetMesh() )
gid_io.FinalizeMesh()
gid_io.InitializeResults( 0.0, model_part.GetMesh() )

for step in range(1,nsteps):
    

    time = Dt*step
    model_part.CloneTimeStep(time)
    print "line49"
    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        

        solver.Solve()
                #print the results
        gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes, time, 0 )
        
        #gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time)
    

gid_io.FinalizeResults()
          
        

