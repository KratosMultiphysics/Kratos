# Import libraries
import numpy
import scipy
import scipy.sparse
import scipy.sparse.linalg
import time as timemodule                   # Import time library as timemodule (avoid interferences with "time" var)
import json                                 # Encoding library (for data exchange)
# Import utilities
import connectivity_mapper                  # Auxiliary matching meshes communicator
import residual_definitions                 # Residual definitions
import mvqn_strategy                        # MultiVector Quasi-Newton method strategy
import jfnk_strategy                        # Jacobian Free Newton-Krylov method strategy
import relaxation_strategy                  # Relaxation strategy
import KratosMultiphysics                   # Kratos core
import KratosMultiphysics.FSIApplication    # FSIApplication (contains the 3D mapper)
# Import solvers
import FluidProblemClass                    # Fluid problem solver
import SolidProblemClass                    # Solid problem solver
# Import ProjectParameters
import ProjectParametersFluid               # Fluid ProjectParameters file
import ProjectParametersSolid               # Solid ProjectParameters file  


############################################################################################################

# Initial tests
if ProjectParametersFluid.domain_size != ProjectParametersSolid.domain_size:
    raise("ERROR: Different working dimensions among subdomains!")
if ProjectParametersFluid.Dt != ProjectParametersSolid.time_step:
    raise("ERROR: Different time step among subdomains!")
if ProjectParametersFluid.nsteps != ProjectParametersSolid.nsteps:
    raise("ERROR: Different number of time steps among subdomains!")
if ProjectParametersFluid.max_time != ProjectParametersSolid.end_time:
    raise("ERROR: Different final time among subdomains!")
if ProjectParametersFluid.output_time != ProjectParametersSolid.GiDWriteFrequency:
    raise("ERROR: Different output time among subdomains!")

# Stepping and time settings
domain_size = ProjectParametersFluid.domain_size
Dt = ProjectParametersFluid.Dt
Nsteps = ProjectParametersFluid.nsteps
final_time = ProjectParametersFluid.max_time
output_time = ProjectParametersFluid.output_time
time = ProjectParametersFluid.Start_time

# Solid and fluid problem construction
D1_problem = FluidProblemClass.FluidProblem(ProjectParametersFluid)  ## D1_problem corresponds to the fluid domain
D2_problem = SolidProblemClass.SolidProblem(ProjectParametersSolid)  ## D2_problem corresponds to the solid domain

# Solid and fluid problem initialization
print("Fluid and solid problems initialization...")
D1_problem.Initialize()
D2_problem.Initialize()
print("Fluid and solid problems initialization finished.")

#~ # Add axiliar FSI mapper variables
#~ import NonConformant_OneSideMap
#~ NonConformant_OneSideMap.AddVariables(D1_problem.fluid_model_part, D2_problem.model_part)

# Conditions initialization
print("Fluid and solid problem dependent conditions set up starts...")
D1_problem.InitializeConditions()
D2_problem.InitializeConditions()
print("Fluid and solid problem dependent conditions finished.")

# Fluid interface is taken as reference interface
Interface_pb_size=D1_problem.interface_nodes

# Interface communicator construction
if domain_size == 2:
    # Fetch the solid and fluid interface nodes
    fluid_interface_nodes=D1_problem.interface_nodes_vec
    solid_interface_nodes=D2_problem.interface_nodes_vec
    
    # Check whether the interface nodes match or not
    if D1_problem.interface_nodes != D2_problem.interface_nodes:
        raise("ERROR: Different number of interface nodes among subdomains!. Non-conformant 2D mapper not implemented yet.")
        
    # Construct the 2D conformant interface mapper
    print("2D interface communicator construction starts...")
    wet_interface_comm = connectivity_mapper.interface_communicator(fluid_interface_nodes,solid_interface_nodes)
    print("2D interface communicator successfully constructed.")
    
#~ elif domain_size == 3:
    #~ print("3D non-conformant interface mapper construction starts...")
    #~ wet_interface_comm = NonConformant_OneSideMap.NonConformant_OneSideMap(D1_problem.fluid_model_part,D2_problem.model_part, 1.0, 15)
    #~ print("3D non-conformant interface successfully constructed.")
    # TAL COMO ESTÁ EL RESIDUO SÓLO VA A FUNCIONAR PARA 2D. LAS LLAMADAS DEL MAPPER SE LLAMAN IGUAL, PERO LOS ARGUMENTOS NO SON LOS MISMOS.

# Output initialization
#### VERY IMPORTANT ####
# Note that the solid problem is firstly initialized. Otherwise, the fluid output files dissapear as soon as the solid output is initialized.
# TODO: Fix this issue
print("Output initialization...")
D2_problem.InitializeOutput()
D1_problem.InitializeOutput()
print("Output initialization finished.")

# Interface residual construction
coupling_algorithm = "NeumannNeumann"

print("Interface residual construction starts...")
if coupling_algorithm == "DirichletNeumann":
    residual = residual_definitions.DirichletNeumannResidual(D1_problem,D2_problem,wet_interface_comm)
    print("Dirichlet-Neumann residual constructed.")
    
elif coupling_algorithm == "NeumannNeumann":
    residual = residual_definitions.NeumannNeumannResidual(D1_problem,D2_problem,wet_interface_comm)
    print("Neumann-Neumann residual constructed.")
    
# Interface strategy construction
coupling_strategy = "JFNK"

print("Interface coupling strategy construction starts...")
if coupling_strategy == "MVQN":
    interface_strategy = mvqn_strategy.MultiVectorQuasiNewtonStrategy(Interface_pb_size,domain_size,D1_problem,D2_problem,wet_interface_comm,residual)
    print("MultiVector Quasi-Newton strategy constructed.")
elif coupling_strategy == "JFNK":
    interface_strategy = jfnk_strategy.JacobianFreeNewtonKrylovStrategy(Interface_pb_size,domain_size,D1_problem,D2_problem,wet_interface_comm,residual)
    print("Jacobian Free Newton-Krylov strategy constructed.")
elif coupling_strategy == "Relaxation":
    if coupling_algorithm == "DirichletNeumann":
        interface_strategy = relaxation_strategy.RelaxationStrategy(Interface_pb_size,domain_size,D1_problem,D2_problem,wet_interface_comm,residual)
        print("Relaxation strategy constructed.")
    elif coupling_algorithm == "NeumannNeumann":
        raise("ERROR: Relaxation strategy must be used with Dirichlet-Neumann algorithm.")
else:
    raise("ERROR: Interface strategy not implemented yet!")

# Output files 
filename = "MAIN_FILE_FSI_"+coupling_algorithm+"_"+coupling_strategy+".log"

# .log File creation to store the iterations evolution
with open(filename, 'w') as file:    # w: only writting (an existing file w/ same name would be erased)
    file.write("Interface problem size: "+str(Interface_pb_size)+"\n")
    file.write("Interface residual size: "+str(Interface_pb_size*domain_size)+"\n"+"\n")
    file.close()

# NL solver parameters
max_nl_iterations = 50
nl_tol = 1e-5
    
iteration_guess_value = 0.0001*numpy.ones(Interface_pb_size*domain_size,dtype='float')       # Interface unknown guess (it might be velocity or fluxes depending on the type of coupling)                            

out = 0
step = 0

print("COUPLED PROBLEM RESOLUTION STARTS...")
print("Interface problem size: ",Interface_pb_size)

while(time <= final_time):
    
    time = time + Dt
    step = step + 1
    convergence = False
    
    D1_problem.ExecuteInitializeSolutionStep(time)
    D2_problem.ExecuteInitializeSolutionStep(time,step)
    
    interface_strategy.ExecuteInitializeSolutionStep()
   
    print("STEP = ", step)
    print("TIME = ", time)
    
    with open(filename, 'a') as file:
        file.write("STEP: "+str(step)+"\n")
        file.write("TIME: "+str(time)+"\n")
        file.close()

    for nl_it in range(1,max_nl_iterations+1): # Problem linearization solving several GMRES minimization problems.
        print("    NL-ITERATION ",nl_it,"STARTS.")
        
        print("     Residual computation starts...")
        vel_residual = residual.ComputeResidual(iteration_guess_value)
        nl_res_norm = scipy.linalg.norm(vel_residual)
        print("     Residual computation finished. |res|=",nl_res_norm)
        
        ### CONVERGENCE ACHIEVED ###
        if nl_res_norm < nl_tol: 
            
            convergence = True
            
            print("    CONVERGENCE ACHIEVED")
            print("    Total non-linear iterations: ",nl_it," NL residual norm: ",nl_res_norm)
            
            with open(filename, 'a') as file:
                file.write("    Non-linear iteration summary:\n")
                file.write("        nl_it: "+str(nl_it)+"\n")
                file.write("        nl_res_norm: "+str(nl_res_norm)+"\n"+"\n")
                file.write("    CONVERGENCE ACHIEVED\n"+"\n")
                file.close()      
            
            break
        
        ### CONVERGENCE NOT ACHIEVED ###
        else:
            
            iteration_corrected_value = interface_strategy.InterfaceSolutionUpdate(step,nl_it,iteration_guess_value,vel_residual)
            iteration_guess_value = numpy.copy(iteration_corrected_value)
        
            with open(filename, 'a') as file:
                file.write("    Non-linear iteration summary:\n")
                file.write("        nl_it: "+str(nl_it)+"\n")
                file.write("        nl_res_norm: "+str(nl_res_norm)+"\n"+"\n")
                file.close()
                    
    # Solve the mesh movement
    solid_interface_disp = D2_problem.GetInterfaceDisplacement()
    solid_interface_disp_comm = wet_interface_comm.StructureToFluid_VectorMap(solid_interface_disp)
    
    D1_problem.SolveMesh(solid_interface_disp_comm)
    
    # Print results
    if(output_time <= out):
        out = 0

    out = out + Dt
    
    D1_problem.ExecuteFinalizeSolutionStep(time,output_time,out)
    D2_problem.ExecuteFinalizeSolutionStep()
    
    interface_strategy.ExecuteFinalizeSolutionStep()
        
D1_problem.ExecuteFinalize()
D2_problem.ExecuteFinalize()
    
print("COUPLED PROBLEM SOLVED.")
