##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2
import mpi;

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs' ##kratos_root/libs
#kratos_libs_path = 'C:/kratosR1/libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
import sys

sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

print "importing KratosMetisApplication ..."
sys.path.append(kratos_applications_path + '/metis_application/python_scripts')
from KratosMetisApplication import *
metis_application = KratosMetisApplication()
kernel.AddApplication(metis_application)
print "KratosMetisApplication sucessfully imported"

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

kernel.InitializeApplication(metis_application)


## from now on the order is not anymore crucial
##################################################################
##################################################################


#defining a model part
model_part = ModelPart("FluidPart");
model_part.AddNodalSolutionStepVariable(VELOCITY)
model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
use_multi_file = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
 
gid_io = GidIO("test",gid_mode,use_multi_file,deformed_mesh_flag, write_conditions)

#print kernel

print "line 58"

number_of_partitions = mpi.size #we set it equal to the number of processors       
#metis_partitioning_process = MetisPartitioningProcess(model_part, gid_io, number_of_partitions)
contact_indices = IndicesVector()
contact_indices[:] = [1,2,3,6,7,11,12]
metis_partitioning_process = MetisContactPartitioningProcess(model_part, gid_io, number_of_partitions, contact_indices)

print "line 62"
print GetRank(), ": i am : ... "

import time
time.sleep(2.0)

metis_partitioning_process.Execute()

print "line 65"
print GetRank(), ": metis_partitioning_process finised"

print model_part



gid_io.InitializeMesh( GetRank()  )
gid_io.WriteMesh( model_part.GetMesh() )
gid_io.FinalizeMesh()
gid_io.InitializeResults( GetRank(), model_part.GetMesh() )
gid_io.WriteNodalResults( PARTITION_INDEX, model_part.Nodes, 0, 0 )
gid_io.FinalizeResults();  

##################################################################
## command line:
## mpirun -np 2 /usr/bin/mpipython square_domain.py
## where 2 is the mumber of threads
##################################################################


