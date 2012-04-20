##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

#defining a model part
model_part = ModelPart("FluidPart");  
model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
model_part.AddNodalSolutionStepVariable(TEMPERATURE)

#adding of Variables to Model Part should be here
import petsc_structural_solver_static
petsc_structural_solver_static.AddVariables(model_part)

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary    # or GiDPostMode.GiD_PostAscii
use_multi_file = MultiFileFlag.MultipleFiles    # or MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed    # or WriteDeformedMeshFlag.WriteUndeformedMesh
write_conditions = WriteConditionsFlag.WriteElementsOnly   # or WriteConditionsFlag.WriteConditions

gid_io = GidIO("cantilever2d",gid_mode,use_multi_file,deformed_mesh_flag, write_conditions)
#gid_io.ReadModelPart(model_part)

metis_partitioning_process = MetisPartitioningProcess(model_part, gid_io)

if  GetRank() == 0 :
    print gid_io


metis_partitioning_process.Execute()

print model_part


mesh_name = GetRank();
print "mesh_name", mesh_name
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()

print model_part


#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
petsc_structural_solver_static.AddDofs(model_part)

#creating a fluid solver object
solver = petsc_structural_solver_static.StaticStructuralSolver(model_part,domain_size)
#pILUPrecond = ILU0Preconditioner() 
#solver.structure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
#solver.structure_linear_solver =  PetscSolver()
solver.Initialize()
solver.solver.space_utils = PetscSparseSpace()
solver.solver.scheme = PetscResidualBasedIncrementalUpdateStaticScheme()
solver.solver.builder_and_solver = PetscResidualBasedEliminationBuilderAndSolver()

model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
print "Linear elastic model selected"


(solver).SetEchoLevel(2);

Dt = 0.001
nsteps = 10

gid_io.InitializeResults( mesh_name , model_part.GetMesh() )

for step in range(0,nsteps):

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 2):
        solver.Solve()

#        for node in model_part.Nodes :
#            print "rank : ", mesh_name, "node", node.Id, "equation_id", node.GetSolutionStepValue(TEMPERATURE)


    #print the results
    #gid_io.InitializeResults(mesh_name,(model_part).GetMesh())
gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,step,0)
gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,step,0)
gid_io.WriteNodalResults(REACTION,model_part.Nodes,step,0)
gid_io.WriteNodalResults(POSITIVE_FACE_PRESSURE,model_part.Nodes,step,0)
gid_io.Flush()
    #gid_io.FinalizeResults()
    #gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time,domain_size)
print "finito"

gid_io.FinalizeResults();  

          
##################################################################
## command line:
## mpirun.mpich -np 2 /usr/bin/mpipython cantilever2dstatic.py
## where 2 is the mumber of threads
##################################################################

        

