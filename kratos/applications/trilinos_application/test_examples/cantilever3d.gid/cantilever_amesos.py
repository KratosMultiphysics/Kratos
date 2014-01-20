from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
#
#
# setting the domain size for the problem to be solved
domain_size = 3

import sys
kratos_path = '../../../..'
kratos_benchmarking_path = '../../../../benchmarking'
sys.path.append(kratos_path)
sys.path.append(kratos_benchmarking_path)

from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *

print("i am ", mpi.rank, " of ", mpi.size)

#
#

# defining a model part
model_part = ModelPart("FluidPart")

# adding of Variables to Model Part should be here when the "very fix container will be ready"
import trilinos_structural_solver_static
trilinos_structural_solver_static.AddVariables(model_part)
model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

# reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("cantilever3d", gid_mode, multifile, deformed_mesh_flag, write_conditions)
# gid_io.ReadModelPart(model_part)

number_of_partitions = mpi.size  # we set it equal to the number of processors
print("number_of_partitions", number_of_partitions)
partitioner = MetisPartitioningProcess(model_part, gid_io, number_of_partitions, domain_size)
partitioner.Execute()
print("GetRank()", GetRank())

mesh_name = mpi.rank
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh((model_part).GetMesh())
gid_io.FinalizeMesh()


print("mesh_name =", mesh_name)
print(model_part)
print(model_part.Properties)

# writing the mesh
# gid_io.WriteUndeformedMesh(model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);

# the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

# importing the solver files
trilinos_structural_solver_static.AddDofs(model_part)

# creating a fluid solver object
solver = trilinos_structural_solver_static.StaticStructuralSolver(model_part, domain_size)
# pILUPrecond = ILU0Preconditioner()
# solver.structure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
solver_parameters = ParameterList()
# solver_parameters.set("PrintTiming", True);
# solver_parameters.set("PrintStatus", True);
solver.structure_linear_solver = AmesosSolver("Superludist", solver_parameters)
model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D())
print("Linear elastic model selected")

solver.Initialize()
(solver.solver).SetEchoLevel(2)

Dt = 0.001
nsteps = 10

gid_io.InitializeResults(mesh_name, (model_part).GetMesh())

for step in range(0, nsteps):
    print("line49")

    time = Dt * step
    model_part.CloneTimeStep(time)

    print(time)
    # print model_part.ProcessInfo()[TIME]

    # solving the fluid problem
    if(step > 3):
        solver.Solve()

#    if(step > 4):

        # print the results
        gid_io.WriteNodalResults(DISPLACEMENT, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(REACTION, model_part.Nodes, time, 0)
#        gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time,domain_size)

gid_io.FinalizeResults()
