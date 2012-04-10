##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################

kratos_path = '../../../..'
import sys
sys.path.append(kratos_path)

from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *

##################################################################

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here
import trilinos_structural_solver_static
trilinos_structural_solver_static.AddVariables(model_part)
model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("cantilever2d",gid_mode,multifile,deformed_mesh_flag, write_conditions)


number_of_partitions = mpi.size #we set it equal to the number of processors
print "number_of_partitions", number_of_partitions
partitioner = MetisPartitioningProcess(model_part, gid_io, number_of_partitions, domain_size);
partitioner.Execute()

mesh_name = mpi.rank
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
trilinos_structural_solver_static.AddDofs(model_part)

#creating a fluid solver object
solver = trilinos_structural_solver_static.StaticStructuralSolver(model_part,domain_size)
##pILUPrecond = ILU0Preconditioner() 
##solver.structure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)

model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
print "Linear elastic model selected"


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
gid_io.FinalizeResults()
print "finito"


          
        

