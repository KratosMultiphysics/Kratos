from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import mpi
import fluid_only_var


#
#
# setting the domain size for the problem to be solved
domain_size = fluid_only_var.domain_size


kratos_benchmarking_path = fluid_only_var.kratos_path + \
    '/benchmarking'  # kratos_root/benchmarking
import sys
sys.path.append(kratos_benchmarking_path)
import benchmarking

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *


# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")

#
# importing the solvers needed
SolverType = fluid_only_var.SolverType
if(SolverType == "FractionalStep"):
    import fractional_step_solver
    fractional_step_solver.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import trilinos_monolithic_solver_eulerian
    trilinos_monolithic_solver_eulerian.AddVariables(fluid_model_part)
else:
    raise "solver type not supported: options are FractionalStep - Monolithic"

# introducing input file name
input_file_name = fluid_only_var.problem_name

# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(
    input_file_name,
    gid_mode,
    multifile,
    deformed_mesh_flag,
    write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
# model_part_io_fluid.ReadModelPart(fluid_model_part)


mpi.world.barrier
print("AAAPAPAPAPAPAAPAPAPAPAPAPAPAPAPAPAPAP")

number_of_partitions = mpi.size  # we set it equal to the number of processors
print("number_of_partitions", number_of_partitions)
partitioner = MetisPartitioningProcess(
    fluid_model_part,
    model_part_io_fluid,
    number_of_partitions,
    domain_size)
partitioner.Execute()
print("GetRank()", GetRank())

# mesh to be printed
mesh_name = mpi.rank
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh((fluid_model_part).GetMesh())
gid_io.FinalizeMesh()
gid_io.Flush()

# mpi.world.barrier()
# prova


# for node in fluid_model_part.Nodes:
# print node.X
# print node.Y
#

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(2)

# adding dofs
if(SolverType == "FractionalStep"):
    fractional_step_solver.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    trilinos_monolithic_solver_eulerian.AddDofs(fluid_model_part)


# select here the laplacian form!!!!!!!!!!!!!!!!!
laplacian_form = fluid_only_var.laplacian_form
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

# check to ensure that no node has zero density or pressure
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print("node ", node.Id, " has zero density!")
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print("node ", node.Id, " has zero viscosity!")
        raise 'node with zero VISCOSITY found'

# creating the solvers
# fluid solver
if(SolverType == "FractionalStep"):
    fluid_solver = fractional_step_solver.IncompressibleFluidSolver(
        fluid_model_part, domain_size)
    fluid_solver.laplacian_form = laplacian_form
    # standard laplacian form
    fluid_solver.predictor_corrector = fluid_only_var.predictor_corrector
    fluid_solver.max_press_its = fluid_only_var.max_press_its
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian"):
    fluid_solver = trilinos_monolithic_solver_eulerian.MonolithicSolver(
        fluid_model_part, domain_size)
    oss_swith = fluid_only_var.use_oss
    dynamic_tau = fluid_only_var.dynamic_tau
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith)
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau)

#
# defining the linear solver
# aztec_parameters = ParameterList()
# aztec_parameters.set("AZ_solver","AZ_gmres");
# aztec_parameters.set("AZ_kspace",100);
# aztec_parameters.set("AZ_output",32);
#
# preconditioner_type = "Amesos"
# preconditioner_parameters = ParameterList()
# preconditioner_parameters.set("amesos: solver type", "Amesos_Klu");
#
# preconditioner_type = "ILU"
# preconditioner_parameters = ParameterList()
#
# overlap_level = 3
# nit_max = 300
# tol = 1e-6
#
# fluid_solver.linear_solver =  AztecSolver(aztec_parameters,preconditioner_type,preconditioner_parameters,tol,nit_max,overlap_level);
#
# fluid_solver.buildertype = "standard"
    fluid_solver.Initialize()


print("fluid solver created")

# settings to be changed
Dt = fluid_only_var.Dt
full_Dt = Dt
initial_Dt = 0.001 * full_Dt  # 0.05 #0.01
final_time = fluid_only_var.max_time
output_step = fluid_only_var.output_step

out = 0


gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())

#
# for node in fluid_model_part.Nodes:
# print node.X
# print node.Y
# aaa

#
# back_node = NodeFinder(fluid_model_part.Nodes , 6.199, 4.0, 0.0)
# print back_node
#
#


time = 0.0
step = 0
while(time < final_time):

    if(step < 5):
        Dt = initial_Dt
    else:
        Dt = full_Dt

    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    if(step >= 3):
        fluid_solver.Solve()

# gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
# gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
# BenchmarkCheck(time, back_node)

    if(out == output_step):
        gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(AIR_PRESSURE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, time, 0)
        # gid_io.WriteNodalResults(DISPLACEMENT,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(
            MESH_VELOCITY,
            fluid_model_part.Nodes,
            time,
            0)
        gid_io.WriteNodalResults(IS_STRUCTURE, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(IS_BOUNDARY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(IS_POROUS, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(
            IS_FREE_SURFACE,
            fluid_model_part.Nodes,
            time,
            0)
        # gid_io.PrintOnGaussPoints(THAWONE,fluid_model_part,time)
        # gid_io.PrintOnGaussPoints(THAWTWO,fluid_model_part,time)
        gid_io.WriteNodalResults(ADVPROJ, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DIVPROJ, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DENSITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DENSITY_AIR, fluid_model_part.Nodes, time, 0)
       # gid_io.WriteNodalResults(NODAL_H,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VISCOSITY, fluid_model_part.Nodes, time, 0)
        gid_io.Flush()

        mpi.world.barrier

        out = 0

    out = out + 1
    step = step + 1

gid_io.FinalizeResults()
