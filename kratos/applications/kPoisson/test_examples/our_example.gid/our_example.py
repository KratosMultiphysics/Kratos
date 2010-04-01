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
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

print "before importing kratos"

#importing Kratos main library
from Kratos import *
print "Kratos library imported"
kernel = Kernel()   #defining kernel
print "kernel created"
#importing applications
import applications_interface
applications_interface.Import_PoissonApplication = True

applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosR1PoissonApplication import *

#defining a model part
model_part = ModelPart("PoissonPart");  

#adding of Variables to Model Part
import static_poisson_solver
static_poisson_solver.AddVariables(model_part)

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("our_example",gid_mode,multifile,deformed_mesh_flag,write_conditions)
gid_io.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#importing the solver files
static_poisson_solver.AddDofs(model_part)
    
#creating a solver object
solver = static_poisson_solver.StaticPoissonSolver(model_part,domain_size)
solver.time_order = 1
solver.linear_solver = SkylineLUFactorizationSolver();
solver.echo_level = 0
solver.Initialize()

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())

solver.Solve()

gid_io.WriteNodalResults(DUMMY_UNKNOWN,model_part.Nodes,0,0)
gid_io.FinalizeResults()
