import fluid_only_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = fluid_only_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
import sys
sys.path.append(kratos_benchmarking_path)
import benchmarking

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *


def BenchmarkCheck(time, model_part):
    max_press = 0.0; 
    min_press = 0.0;
    vel2min = 10000.0;
    id_min_vel = 0
    x_min_vel = 0.0
    y_min_vel = 0.0
    for node in model_part.Nodes:
        press = node.GetSolutionStepValue(PRESSURE);
        if(press > max_press):
            max_press = press
        elif(press < min_press):
             min_press = press

        x = node.X
        y = node.Y
        vel = node.GetSolutionStepValue(VELOCITY)
        vel2 = vel[0]**2 + vel[1]**2
        if(x > 0.1 and x<0.9 and y>0.1 and y<0.9):
            if(vel2 < vel2min):
                vel2min = vel2
                id_min_vel = node.Id
                x_min_vel = node.X
                y_min_vel = node.Y
            
        
    benchmarking.Output(time, "Time")
    benchmarking.Output(min_press, "minimum pressure", 0.00001)
    benchmarking.Output(max_press, "maximum pressure", 0.00001)
    benchmarking.Output(id_min_vel, "Id of the node with minimum velocity norm", 0.0)
    benchmarking.Output(x_min_vel, "coord x minimum velocity norm", 0.0)
    benchmarking.Output(y_min_vel, "coord y minimum velocity norm", 0.0)

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

#########select here the laplacian form!!!!!!!!!!!!!!!!!
laplacian_form = 1

if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

##check to ensure that no node has zero density or pressure
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise 'node with zero VISCOSITY found'

#creating the solvers
#fluid solver
if(SolverType == "fractional_step"):
    fluid_solver = fractional_step_solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
    fluid_solver.laplacian_form = laplacian_form; #standard laplacian form
    fluid_solver.predictor_corrector = False
    fluid_solver.max_press_its = fluid_only_var.max_press_its
    fluid_solver.velocity_linear_solver = SkylineLUFactorizationSolver();
    fluid_solver.pressure_linear_solver = SkylineLUFactorizationSolver();
    
    fluid_solver.Initialize()
elif(SolverType == "pressure_splitting"):
    fluid_solver = decoupled_solver_eulerian.\
		  DecoupledSolver(fluid_model_part,domain_size)
    oss_switch = fluid_only_var.use_oss
    dynamic_tau = fluid_only_var.dynamic_tau
##    pPrecond = ILU0Preconditioner()
    pPrecond = DiagonalPreconditioner()
    fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pPrecond)
##    fluid_solver.linear_solver =  SuperLUSolver()
    fluid_solver.rel_vel_tol = 1e-4
    fluid_solver.abs_vel_tol = 1e-6
    fluid_solver.rel_pres_tol = 1e-4
    fluid_solver.abs_pres_tol = 1e-6
    fluid_solver.use_inexact_newton = False
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch)
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau)
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian"): 
    fluid_solver = monolithic_solver_eulerian.MonolithicSolver(fluid_model_part,domain_size)
    oss_switch = fluid_only_var.use_oss
    dynamic_tau = fluid_only_var.dynamic_tau
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch);				
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian_compressible"): 
    fluid_solver = monolithic_solver_eulerian_compressible.MonolithicSolver(fluid_model_part,domain_size)
    oss_switch = fluid_only_var.use_oss
    dynamic_tau = fluid_only_var.dynamic_tau
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch);				
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.Initialize()


print "fluid solver created"

#settings to be changed
Dt = fluid_only_var.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01
final_time = fluid_only_var.max_time
output_step = fluid_only_var.output_step

out = 0


#mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( fluid_model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())


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
        BenchmarkCheck(time, fluid_model_part)
        

    if(out == output_step):
        if(SolverType == "fractional_step" or SolverType == "pressure_splitting"):
            gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
        else:
            gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(AIR_PRESSURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(WATER_PRESSURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
            #gid_io.WriteNodalResults(DISPLACEMENT,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(MESH_VELOCITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(IS_STRUCTURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(IS_BOUNDARY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(IS_POROUS,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(IS_FREE_SURFACE,fluid_model_part.Nodes,time,0)
            #gid_io.PrintOnGaussPoints(THAWONE,fluid_model_part,time)
            #gid_io.PrintOnGaussPoints(THAWTWO,fluid_model_part,time)
            gid_io.WriteNodalResults(ADVPROJ,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(DIVPROJ,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(DENSITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(DENSITY_AIR,fluid_model_part.Nodes,time,0)
           ## gid_io.WriteNodalResults(NODAL_H,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(VISCOSITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(SOUND_VELOCITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(AIR_SOUND_VELOCITY,fluid_model_part.Nodes,time,0)

        out = 0

    out = out + 1
    step = step + 1
      
gid_io.FinalizeResults()
          
        

