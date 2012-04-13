import problem_settings

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = problem_settings.domain_size

from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  

#############################################
##importing the solvers needed
SolverType = problem_settings.SolverType
if(SolverType == "fractional_step"):
    import fractional_step_solver
    fractional_step_solver.AddVariables(fluid_model_part)
elif(SolverType == "pressure_splitting"):
    import decoupled_solver_eulerian
    decoupled_solver_eulerian.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import monolithic_solver_eulerian
    monolithic_solver_eulerian.AddVariables(fluid_model_part)
    fluid_model_part.AddNodalSolutionStepVariable(YIELD_STRESS)
    fluid_model_part.AddNodalSolutionStepVariable(TAU)
    fluid_model_part.AddNodalSolutionStepVariable(MU)
    fluid_model_part.AddNodalSolutionStepVariable(EQ_STRAIN_RATE)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    import monolithic_solver_eulerian_compressible
    monolithic_solver_eulerian_compressible.AddVariables(fluid_model_part)
else:
    raise Exception("solver type not supported: options are fractional_step - \
	pressure_splitting - monolithic_solver_eulerian - \
	monolithic_solver_eulerian_compressible")

#introducing input file name
input_file_name = problem_settings.problem_name

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile
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
laplacian_form = problem_settings.laplacian_form 
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

##check to ensure that no node has zero density or pressure
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise Exception('node with zero density found')
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise Exception('node with zero VISCOSITY found')

#creating the solvers
#fluid solver
if(SolverType == "fractional_step"):
    fluid_solver = fractional_step_solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
    fluid_solver.laplacian_form = laplacian_form; #standard laplacian form
    fluid_solver.predictor_corrector = problem_settings.predictor_corrector
    fluid_solver.max_press_its = problem_settings.max_press_its
    fluid_solver.Initialize()
elif(SolverType == "pressure_splitting"):
    fluid_solver = decoupled_solver_eulerian.DecoupledSolver(fluid_model_part,domain_size)
##    pPrecond = ILU0Preconditioner()
    pPrecond = DiagonalPreconditioner()
    fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pPrecond)
##    fluid_solver.linear_solver =  SuperLUSolver()
    fluid_solver.rel_vel_tol = 1e-4
    fluid_solver.abs_vel_tol = 1e-6
    fluid_solver.rel_pres_tol = 1e-4
    fluid_solver.abs_pres_tol = 1e-6
    fluid_solver.use_inexact_newton = False
    fluid_solver.dynamic_tau = problem_settings.dynamic_tau
    fluid_solver.oss_switch  = int(problem_settings.use_oss)
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian"): 
    fluid_solver = monolithic_solver_eulerian.MonolithicSolver(fluid_model_part,domain_size)
    fluid_solver.dynamic_tau = problem_settings.dynamic_tau
    fluid_solver.oss_switch  = int(problem_settings.use_oss)
    fluid_solver.regularization_coef = problem_settings.m_coef
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian_compressible"): 
    fluid_solver = monolithic_solver_eulerian_compressible.MonolithicSolver(fluid_model_part,domain_size)
    fluid_solver.dynamic_tau = problem_settings.dynamic_tau
    fluid_solver.oss_switch  = int(problem_settings.use_oss)
    fluid_solver.Initialize()


print "fluid solver created"

#settings to be changed
Dt = problem_settings.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01
final_time = problem_settings.max_time
output_step = problem_settings.output_step

out = 0


#mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( fluid_model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())

# Write .post.list file (GiD postprocess list)
f = open(problem_settings.problem_name+'.post.lst','w')
f.write('Single\n')
f.write(problem_settings.problem_name+'.post.bin\n')
f.close()

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
            gid_io.WriteNodalResults(EXTERNAL_PRESSURE,fluid_model_part.Nodes,time,0)
##            gid_io.PrintOnGaussPoints(TEMPERATURE,fluid_model_part,time)
##            gid_io.PrintOnGaussPoints(AUX_INDEX,fluid_model_part,time)
            gid_io.PrintOnGaussPoints(EQ_STRAIN_RATE,fluid_model_part,time)
	    gid_io.WriteNodalResults(YIELD_STRESS,fluid_model_part.Nodes,time,0)
	    gid_io.PrintOnGaussPoints(MU,fluid_model_part,time)
            gid_io.PrintOnGaussPoints(TAU,fluid_model_part,time)
        out = 0

    out = out + 1
    step = step + 1
      
gid_io.FinalizeResults()
          
        

