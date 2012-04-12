import fluid_only_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = fluid_only_var.domain_size

kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
import sys
sys.path.append(kratos_benchmarking_path)

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *
import benchmarking

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  

#############################################
##importing the solvers needed
SolverType = fluid_only_var.SolverType
if(SolverType == "fractional_step"):
    import fractional_step_solver
    fractional_step_solver.AddVariables(fluid_model_part)
elif(SolverType == "pressure_splitting"):
    import decoupled_solver_eulerian
    decoupled_solver_eulerian.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import monolithic_solver_eulerian
    monolithic_solver_eulerian.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    import monolithic_solver_eulerian_compressible
    monolithic_solver_eulerian_compressible.AddVariables(fluid_model_part)
else:
    raise "solver type not supported: options are fractional_step - \
	pressure_splitting - monolithic_solver_eulerian - \
	monolithic_solver_eulerian_compressible"

#introducing input file name
input_file_name = fluid_only_var.problem_name

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

##adding dofs
if(SolverType == "fractional_step"):
    fractional_step_solver.AddDofs(fluid_model_part)
elif(SolverType == "pressure_splitting"):
    decoupled_solver_eulerian.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    monolithic_solver_eulerian.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    monolithic_solver_eulerian_compressible.AddDofs(fluid_model_part)


#do mesh improvement
reconnector = TetrahedraReconnectUtility()
reconnector.EvaluateQuality(fluid_model_part)

#mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( fluid_model_part.GetMesh() )
gid_io.FinalizeMesh()

time=0.0
gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())
gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
gid_io.FinalizeResults()


          
        

