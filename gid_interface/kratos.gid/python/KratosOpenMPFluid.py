#################################################################
##################################################################
#import the configuration data as read from the GiD
import ProjectParameters



def PrintResults(model_part):
    print "Writing results. Please run Gid for viewing results of analysis."
    for variable_name in ProjectParameters.nodal_results:
        gid_io.WriteNodalResults(variables_dictionary[variable_name],model_part.Nodes,time,0)
    for variable_name in ProjectParameters.gauss_points_results:
        gid_io.PrintOnGaussPoints(variables_dictionary[variable_name],model_part,time)


##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = ProjectParameters.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = ProjectParameters.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_ExternalSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosIncompressibleFluidApplication import *
from KratosExternalSolversApplication import *


## defining variables to be used

variables_dictionary = {"PRESSURE" : PRESSURE,
                        "VELOCITY" : VELOCITY,
                        "REACTION" : REACTION,
                        "DISTANCE" : DISTANCE,}

#defining a model part for the fluid 
fluid_model_part = ModelPart("FluidPart");

if "REACTION" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(REACTION)
if "DISTANCE" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

#############################################
##importing the solvers needed
SolverType = ProjectParameters.SolverType
if(SolverType == "FractionalStep"):
    import incompressible_fluid_solver
    incompressible_fluid_solver.AddVariables(fluid_model_part)
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
    raise "solver type not supported: options are FractionalStep - pressure_splitting - Monolithic"

#introducing input file name
input_file_name = ProjectParameters.problem_name

#reading the fluid part

# initialize GiD  I/O
if ProjectParameters.GiDPostMode == "Binary":
    gid_mode = GiDPostMode.GiD_PostBinary
elif ProjectParameters.GiDPostMode == "Ascii":
    gid_mode = GiDPostMode.GiD_PostAscii
elif ProjectParameters.GiDPostMode == "AsciiZipped":
    gid_mode = GiDPostMode.GiD_PostAsciiZipped
else:
    print "Unknown GiD post mode, assuming Binary"
    gid_mode = GiDPostMode.GiD_PostBinary

if ProjectParameters.GiDWriteMeshFlag == True:
    deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
else:
    deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed

if ProjectParameters.GiDWriteConditionsFlag == True:
    write_conditions = WriteConditionsFlag.WriteConditions
else:
    write_conditions = WriteConditionsFlag.WriteElementsOnly

if ProjectParameters.GiDMultiFileFlag == "Single":
    multifile = MultiFileFlag.SingleFile
elif ProjectParameters.GiDMultiFileFlag == "Multiples":
    multifile = MultiFileFlag.MultipleFiles
else:
    print "Unknown GiD multiple file mode, assuming Single"
    multifile = MultiFileFlag.SingleFile
    
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
if SolverType == "FractionalStep":
    fluid_model_part.SetBufferSize(3)
else:
    fluid_model_part.SetBufferSize(2)

##adding dofs
if(SolverType == "FractionalStep"):
    incompressible_fluid_solver.AddDofs(fluid_model_part)
elif(SolverType == "pressure_splitting"):
    decoupled_solver_eulerian.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    monolithic_solver_eulerian.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    monolithic_solver_eulerian_compressible.AddDofs(fluid_model_part)

#########select here the laplacian form!!!!!!!!!!!!!!!!!
laplacian_form = ProjectParameters.laplacian_form 
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

## Commented out for efficiency, problem type assigns the values automatically
####check to ensure that no node has zero density or pressure
##for node in fluid_model_part.Nodes:
##    if(node.GetSolutionStepValue(DENSITY) == 0.0):
##        print "node ",node.Id," has zero density!"
##        raise 'node with zero density found'
##    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
##        print "node ",node.Id," has zero viscosity!"
##        raise 'node with zero VISCOSITY found'

# decoupled solver schemes: need solver for pressure and solver fol velocity

if ProjectParameters.SolverType in ["pressure_splitting","FractionalStep"]:
    # Velocity preconditioner
    try:
        if ProjectParameters.Velocity_Preconditioner_type == 'Diagonal':
            vel_precond = DiagonalPreconditioner()
        elif ProjectParameters.Velocity_Preconditioner_type == 'ILU0':
            vel_precond = ILU0Preconditioner()
    except AttributeError:
        # If we are using a direct solver Velocity_Preconditioner_type will not
        # exist. This is not an error.
        vel_precond = None

    # Velocity solver
    if ProjectParameters.Velocity_Linear_Solver == "Conjugate gradient":
        if ProjectParameters.Velocity_Preconditioner_type == 'None':
            velocity_linear_solver =  CGSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration)
        else:
            velocity_linear_solver =  CGSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration,vel_precond)
    elif ProjectParameters.Velocity_Linear_Solver == "BiConjugate gradient stabilized":
        if ProjectParameters.Velocity_Preconditioner_type == 'None':
            velocity_linear_solver =  BICGSTABSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration)
        else:
            velocity_linear_solver =  BICGSTABSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration,vel_precond)
    elif ProjectParameters.Velocity_Linear_Solver == "GMRES":
        if ProjectParameters.Velocity_Preconditioner_type == 'None':
            velocity_linear_solver =  GMRESSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration)
        else:
            velocity_linear_solver =  GMRESSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration,vel_precond)
    elif ProjectParameters.Velocity_Linear_Solver == "Skyline LU factorization":
        velocity_linear_solver = SkylineLUFactorizationSolver()
    elif ProjectParameters.Velocity_Linear_Solver == "Super LU":
        velocity_linear_solver = SuperLUSolver()
    elif ProjectParameters.Velocity_Linear_Solver == "Parallel MKL Pardiso":
        velocity_linear_solver = ParallelMKLPardisoSolver()

    # Pressure preconditioner
    try:
        if ProjectParameters.Pressure_Preconditioner_type == 'Diagonal':
            press_precond = DiagonalPreconditioner()
        elif ProjectParameters.Pressure_Preconditioner_type == 'ILU0':
            press_precond = ILU0Preconditioner()
    except AttributeError:
        # If we are using a direct solver Pressure_Preconditioner_type will not
        # exist. This is not an error.
        press_precond = None

    # Pressure solver
    if ProjectParameters.Pressure_Linear_Solver == "Conjugate gradient":
        if ProjectParameters.Pressure_Preconditioner_type == 'None':
            pressure_linear_solver =  CGSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration)
        else:
            pressure_linear_solver =  CGSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration,press_precond)
    elif ProjectParameters.Pressure_Linear_Solver == "BiConjugate gradient stabilized":
        if ProjectParameters.Pressure_Preconditioner_type == 'None':
            pressure_linear_solver =  BICGSTABSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration)
        else:
            pressure_linear_solver =  BICGSTABSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration,press_precond)
    elif ProjectParameters.Pressure_Linear_Solver == "GMRES":
        if ProjectParameters.Pressure_Preconditioner_type == 'None':
            pressure_linear_solver =  GMRESSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration)
        else:
            pressure_linear_solver =  GMRESSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration,press_precond)
    elif ProjectParameters.Pressure_Linear_Solver == "Skyline LU factorization":
        pressure_linear_solver = SkylineLUFactorizationSolver()
    elif ProjectParameters.Pressure_Linear_Solver == "Super LU":
        pressure_linear_solver = SuperLUSolver()
    elif ProjectParameters.Pressure_Linear_Solver == "Parallel MKL Pardiso":
        pressure_linear_solver = ParallelMKLPardisoSolver()
elif "monolithic_solver_eulerian": # single coupled solver
    # preconditioner
    try:
        if ProjectParameters.Monolithic_Preconditioner_type == 'Diagonal':
            precond = DiagonalPreconditioner()
        elif ProjectParameters.Monolithic_Preconditioner_type == 'ILU0':
            precond = ILU0Preconditioner()
    except AttributeError:
        # If we are using a direct solver Pressure_Preconditioner_type will not
        # exist. This is not an error.
        precond = None

    # solver
    if ProjectParameters.Monolithic_Linear_Solver == "Conjugate gradient":
        if ProjectParameters.Monolithic_Preconditioner_type == 'None':
            monolithic_linear_solver =  CGSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration)
        else:
            monolithic_linear_solver =  CGSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration,precond)
    elif ProjectParameters.Monolithic_Linear_Solver == "BiConjugate gradient stabilized":
        if ProjectParameters.Monolithic_Preconditioner_type == 'None':
            monolithic_linear_solver =  BICGSTABSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration)
        else:
            monolithic_linear_solver =  BICGSTABSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration,precond)
    elif ProjectParameters.Monolithic_Linear_Solver == "GMRES":
        if ProjectParameters.Monolithic_Preconditioner_type == 'None':
            monolithic_linear_solver =  GMRESSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration)
        else:
            monolithic_linear_solver =  GMRESSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration,precond)
    elif ProjectParameters.Monolithic_Linear_Solver == "Skyline LU factorization":
        monolithic_linear_solver = SkylineLUFactorizationSolver()
    elif ProjectParameters.Monolithic_Linear_Solver == "Super LU":
        monolithic_linear_solver = SuperLUSolver()
    elif ProjectParameters.Monolithic_Linear_Solver == "Parallel MKL Pardiso":
        monolithic_linear_solver = ParallelMKLPardisoSolver()

dynamic_tau = ProjectParameters.use_dt_in_stabilization
oss_switch = ProjectParameters.use_orthogonal_subscales
#creating the solvers
#fluid solver
if(SolverType == "FractionalStep"):
    fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
    fluid_solver.max_val_its = ProjectParameters.max_vel_its
    fluid_solver.max_press_its = ProjectParameters.max_press_its
    fluid_solver.laplacian_form = laplacian_form; #standard laplacian form
    fluid_solver.predictor_corrector = ProjectParameters.predictor_corrector
    fluid_solver.use_dt_in_stabilization = False
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch);               
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.vel_toll = ProjectParameters.velocity_relative_tolerance
    fluid_solver.press_toll = ProjectParameters.pressure_relative_tolerance
    fluid_solver.CalculateReactions = ProjectParameters.Calculate_reactions
    # Solver definition
    fluid_solver.velocity_linear_solver = velocity_linear_solver
    fluid_solver.pressure_linear_solver = pressure_linear_solver
    fluid_solver.Initialize()
if(SolverType == "pressure_splitting"):
    fluid_solver = decoupled_solver_eulerian.DecoupledSolver(fluid_model_part,domain_size)
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch);               
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.rel_vel_tol = ProjectParameters.velocity_relative_tolerance
    fluid_solver.abs_vel_tol = ProjectParameters.velocity_absolute_tolerance
    fluid_solver.rel_pres_tol = ProjectParameters.pressure_relative_tolerance
    fluid_solver.abs_pres_tol = ProjectParameters.pressure_absolute_tolerance
    fluid_solver.use_inexact_newton = False
    fluid_solver.max_iter = ProjectParameters.max_iterations
    fluid_solver.CalculateReactions = ProjectParameters.Calculate_reactions
    # Solver definition
    fluid_solver.velocity_linear_solver = velocity_linear_solver
    fluid_solver.pressure_linear_solver = pressure_linear_solver
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian"): 
    fluid_solver = monolithic_solver_eulerian.MonolithicSolver(fluid_model_part,domain_size)
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch);               
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.rel_vel_tol = ProjectParameters.velocity_relative_tolerance
    fluid_solver.abs_vel_tol = ProjectParameters.velocity_absolute_tolerance
    fluid_solver.rel_pres_tol = ProjectParameters.pressure_relative_tolerance
    fluid_solver.abs_pres_tol = ProjectParameters.pressure_absolute_tolerance
    fluid_solver.max_iter = ProjectParameters.max_iterations
    fluid_solver.CalculateReactions = ProjectParameters.Calculate_reactions
    # Solver definition
    fluid_solver.linear_solver = monolithic_linear_solver
    # fluid_solver.pressure_linear_solver = pressure_linear_solver
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian_compressible"): 
    fluid_solver = monolithic_solver_eulerian_compressible.MonolithicSolver(fluid_model_part,domain_size)
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch);               
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.rel_vel_tol = ProjectParameters.velocity_relative_tolerance
    fluid_solver.abs_vel_tol = ProjectParameters.velocity_absolute_tolerance
    fluid_solver.rel_pres_tol = ProjectParameters.pressure_relative_tolerance
    fluid_solver.abs_pres_tol = ProjectParameters.pressure_absolute_tolerance
    fluid_solver.max_iter = ProjectParameters.max_iterations
    fluid_solver.CalculateReactions = ProjectParameters.Calculate_reactions
    # Solver definition
    fluid_solver.linear_solver = monolithic_linear_solver
    # fluid_solver.velocity_linear_solver = velocity_linear_solver
    # fluid_solver.pressure_linear_solver = pressure_linear_solver
    fluid_solver.Initialize()


print "fluid solver created"

#settings to be changed
Dt = ProjectParameters.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01
Nsteps  = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

out = 0


###mesh to be printed (single mesh case)
if ProjectParameters.GiDMultiFileFlag == "Single":
    mesh_name = 0.0
    gid_io.InitializeMesh( mesh_name )
    gid_io.WriteMesh( fluid_model_part.GetMesh() )
    gid_io.FinalizeMesh()

    gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())
    Multifile = False

    # Write .post.list file (GiD postprocess list)
    f = open(ProjectParameters.problem_name+'.post.lst','w')
    f.write('Single\n')
    if ProjectParameters.GiDPostMode == "Binary":
        f.write(ProjectParameters.problem_name+'.post.bin\n')
    elif ProjectParameters.GiDPostMode == "Ascii":
        f.write(ProjectParameters.problem_name+'_0.post.msh\n')
    f.close()
    
else: #  ProjectParameters.GiDMultiFileFlag == "Multiples":
    Multifile = True

    # Initialize .post.list file (GiD postprocess list)
    f = open(ProjectParameters.problem_name+'.post.lst','w')
    f.write('Multiple\n')

Dt      = ProjectParameters.Dt
MaxTime = ProjectParameters.max_time
Nsteps  = ProjectParameters.nsteps

time = 0.0
step = 0
while(time < final_time):

    if(step < 5):
        Dt = initial_Dt
    else:
        Dt = full_Dt
        
    time = time + Dt
    step = step + 1
    fluid_model_part.CloneTimeStep(time)

    print "STEP = ", step
    print "TIME = ", time

    if(step >= 3):
        fluid_solver.Solve()

    if(output_time <= out):
        if Multifile:
            gid_io.InitializeMesh( time )
            gid_io.WriteMesh( fluid_model_part.GetMesh() )
            gid_io.FinalizeMesh()

            gid_io.InitializeResults(time,(fluid_model_part).GetMesh())
        
        PrintResults(fluid_model_part)
        out = 0

        if Multifile:
            gid_io.FinalizeResults()
            if ProjectParameters.GiDPostMode == "Binary":
                f.write(ProjectParameters.problem_name+'_'+str(time)+'.post.bin\n')
            elif ProjectParameters.GiDPostMode == "Ascii":
                f.write(ProjectParameters.problem_name+'_'+str(time)+'.post.msh\n')

    out = out + Dt

if Multifile == False:
    gid_io.FinalizeResults()
    f.close()
          
        
