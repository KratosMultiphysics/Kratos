import mpi
import Kratos_Structural_Application_var

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
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_KratosTrilinosApplication = True
applications_interface.Import_KratosMetisApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosStructuralApplication import *
from KratosIncompressibleFluidApplication import *
from KratosTrilinosApplication import *
from KratosMetisApplication import *
import benchmarking


#defining a model part
model_part = ModelPart("FluidPart");  
input_file_name = Kratos_Structural_Application_var.problem_name

#providing the variable list to the model part
model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

#reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name, gid_mode_flag, use_multifile, deformed_print_flag, write_conditions)



model_part_io_fluid = ModelPartIO(input_file_name)

print "before performing the division"
number_of_partitions = mpi.size #we set it equal to the number of processors
if mpi.rank == 0 :
    partitioner = MetisDivideInputToPartitionsProcess(model_part_io_fluid, number_of_partitions, domain_size);
    partitioner.Execute()

"print division performed"

mpi.world.barrier()

MPICommSetup = SetMPICommunicatorProcess(model_part)
MPICommSetup.Execute()

my_input_filename = input_file_name + "_" + str(mpi.rank)
model_part_io_fluid = ModelPartIO(my_input_filename)
model_part_io_fluid.ReadModelPart(model_part)


print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#adding dofs
for node in model_part.Nodes:
    node.AddDof(DISPLACEMENT_X)
    node.AddDof(DISPLACEMENT_Y)
    node.AddDof(DISPLACEMENT_Z)

model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )



ccc= ParallelFillCommunicator(model_part)
ccc.Execute()

refinement_steps = 1
for i in range(0,refinement_steps):
    print "*******************************************************************"
    for elem in model_part.Elements:
	elem.SetValue(SPLIT_ELEMENT,True)
    
    print "doing refinement step = ",i
    Comm = CreateCommunicator()
    mesh_utility = TrilinosRefineMesh(model_part,Comm)
    refine_on_reference = False
    interpolate_internal_variables = False
    mesh_utility.Local_Refine_Mesh(refine_on_reference,interpolate_internal_variables,domain_size)


    print "-----------------------------------------------------------------"
    ccc= ParallelFillCommunicator(model_part)
    ccc.Execute()
    ##ccc.PrintDebugInfo()

    mpi.world.barrier()

    print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

mesh_name = mpi.rank
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
gid_io.Flush()

gid_io.InitializeResults(mesh_name , model_part.GetMesh())
    
gid_io.FinalizeResults()

          
        

