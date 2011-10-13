#importing MPI ... for this boost 1.35 or superior is needed
try:
 import boost.mpi as mpi
except ImportError:
 import mpi
print "i am ",mpi.rank , " of ",mpi.size

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3

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
applications_interface.Import_StructuralApplication = True
applications_interface.Import_KratosTrilinosApplication = True
applications_interface.Import_KratosMetisApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)
from KratosStructuralApplication import *
from KratosTrilinosApplication import *
from KratosMetisApplication import *

## from now on the order is not anymore crucial
##################################################################
##################################################################

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here when the "very fix container will be ready"
import trilinos_structural_solver_static
trilinos_structural_solver_static.AddVariables(model_part)
model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("cantilever3d",gid_mode,multifile,deformed_mesh_flag, write_conditions)
##gid_io.ReadModelPart(model_part)

number_of_partitions = mpi.size #we set it equal to the number of processors
print "number_of_partitions", number_of_partitions
partitioner = MetisPartitioningProcess(model_part, gid_io, number_of_partitions, domain_size);
partitioner.Execute()
print "GetRank()",GetRank()

mesh_name = mpi.rank
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()


print "mesh_name =", mesh_name
print model_part
print model_part.Properties

#writing the mesh
#gid_io.WriteUndeformedMesh(model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
trilinos_structural_solver_static.AddDofs(model_part)

#creating a fluid solver object
solver = trilinos_structural_solver_static.StaticStructuralSolver(model_part,domain_size)

#defining the linear solver
aztec_parameters = ParameterList()
aztec_parameters.set("AZ_solver","AZ_gmres");
aztec_parameters.set("AZ_kspace",100);
aztec_parameters.set("AZ_output",32);

preconditioner_type = "Amesos"
preconditioner_parameters = ParameterList()
preconditioner_parameters.set("amesos: solver type", "Amesos_Klu");

##preconditioner_type = "ILU"
##preconditioner_parameters = ParameterList()

overlap_level = 3
nit_max = 300
tol = 1e-6

solver.structure_linear_solver =  AztecSolver(aztec_parameters,preconditioner_type,preconditioner_parameters,tol,nit_max,overlap_level);
solver.structure_linear_solver.SetScalingType(AztecScalingType.LeftScaling)
model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
print "Linear elastic model selected"

solver.Initialize()
(solver.solver).SetEchoLevel(2);

Dt = 0.001
nsteps = 10

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())

for step in range(0,nsteps):
    print "line49"

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        solver.Solve()

#    if(step > 4):

        #print the results
        gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
#        gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time,domain_size)

gid_io.FinalizeResults()

          
        

