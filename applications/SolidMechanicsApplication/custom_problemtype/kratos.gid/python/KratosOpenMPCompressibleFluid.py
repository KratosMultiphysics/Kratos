#################################################################
##################################################################
#import the configuration data as read from the GiD
import ProjectParameters

def PrintResults(model_part):
        print "Writing results. Please run Gid for viewing results of analysis."
        for variable_name in ProjectParameters.nodal_results:
            gid_io.WriteNodalResults(varibles_dictionary[variable_name],model_part.Nodes,time,0)
            
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
applications_interface.Import_ConvectionDiffusionApplication = True
applications_interface.Import_PFEMApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosIncompressibleFluidApplication import *
from KratosConvectionDiffusionApplication import *
from KratosExternalSolversApplication import *
from KratosPFEMApplication import *


## defining variables to be used

varibles_dictionary = {"PRESSURE" : PRESSURE,
                       "VELOCITY" : VELOCITY,
                       "TEMPERATURE" : TEMPERATURE,
                       "DENSITY" : DENSITY}

#defining a model part for the fluid 
fluid_model_part = ModelPart("FluidPart");
temperature_model_part = ModelPart("TemperaturePart");

#############################################
##importing the solvers needed
SolverType = ProjectParameters.SolverType
if(SolverType == "FractionalStep_RK4-BE"):
    import runge_kutta_frac_step_comp_solver
    runge_kutta_frac_step_comp_solver.AddVariables(fluid_model_part)
else:
    raise "solver type not supported: only Fractional Step semi-explicit solver available"

#introducing input file name
input_file_name = ProjectParameters.problem_name

import nonlinear_convection_diffusion_solver
nonlinear_convection_diffusion_solver.AddVariables(fluid_model_part) #the nodes are the same 

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
runge_kutta_frac_step_comp_solver.AddDofs(fluid_model_part)
nonlinear_convection_diffusion_solver.AddDofs(fluid_model_part)

NISTTools = NistUtils()

######### check here if the pressure is fixed at least at one point!!!!!!!!!!!!!!!!!
for node in fluid_model_part.Nodes:        
        fixed+=node.IsFixed(PRESSURE)

if (fixed==0):
        print "YOU MUST FIX PRESSURE AT LEAST AT ONE NODE!"
        raise 'Pressure is not fixed anywhere!'

##check to ensure that no node has zero density or density
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise 'node with zero VISCOSITY found'
    if(node.GetSolutionStepValue(CONDUCTIVITY) == 0.0):
        print "node ",node.Id," has zero conductivity!"
        raise 'node with zero CONDUCTIVITY found'
    if(node.GetSolutionStepValue(SPECIFIC_HEAT) == 0.0):
        print "node ",node.Id," has zero specific heat!"
        raise 'node with zero SPECIFIC_HEAT found'

#creating the solvers
#fluid solver
runge_kutta_frac_step_comp_solver.AddDofs(fluid_model_part)
nonlinear_convection_diffusion_solver.AddDofs(fluid_model_part)

fluid_solver = runge_kutta_frac_step_comp_solver.RungeKuttaFracStepCompSolver(model_part,domain_size)
fluid_solver.ReformDofAtEachIteration = False
fluid_solver.Initialize()


temperature_solver = nonlinear_convection_diffusion_solver.ConvectionDiffusionSolver(temperature_model_part,domain_size)
temperature_solver.time_order = 1
    
temperature_solver.ReformDofAtEachIteration = False
temperature_solver.echo_level = 1
temperature_solver.Initialize()



#Linear solver for pressure
if(ProjectParameters.PressureLinearSolver == "SkylineLUFactorization"):
    fluid_solver.linear_solver  =  SkylineLUFactorizationSolver()
elif(ProjectParameters.PressureLinearSolver == "SuperLUSolver"):
    fluid_solver.linear_solver  =  SuperLUSolver()
elif(ProjectParameters.PressureLinearSolver == "CGSolver"):
    pDiagPrecond = DiagonalPreconditioner()
    LST  = ProjectParameters.Linear_Solver_Tolerance
    LSMI = ProjectParameters.Linear_Solver_Max_Iteration  
    fluid_solver.linear_solver  =  CGSolver(LST,LSMI,pDiagPrecond)

print "fluid solver created"
    
#Linear solver for convection-diffusion
if(ProjectParameters.TemperatureLinearSolver == "SkylineLUFactorization"):
    temperature_solver.linear_solver  =  SkylineLUFactorizationSolver()
elif(ProjectParameters.TemperatureLinearSolver == "SuperLUSolver"):
    temperature_solver.linear_solver  =  SuperLUSolver()
elif(ProjectParameters.TemperatureLinearSolver == "CGSolver"):
    pDiagPrecond = DiagonalPreconditioner()
    LST  = ProjectParameters.Linear_Solver_Tolerance
    LSMI = ProjectParameters.Linear_Solver_Max_Iteration  
    temperature_solver.linear_solver  =  CGSolver(LST,LSMI,pDiagPrecond)

print "temperature solver created"    

#settings to be changed
Dt = ProjectParameters.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01
final_time = ProjectParameters.max_time
output_step = ProjectParameters.output_step

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
        Dt=(CFL_time_estimate_process).EstimateTime(CFL, full_Dt)
        print "CFL gave this time step", Dt        
        
    time = time + Dt
    fluid_model_part.CloneTimeStep(time)
    temperature_model_part.CloneTimeStep(time)

    if(step == 3):
        NISTTools.GenerateModelPart(fluid_model_part,temperature_model_part,domain_size);
            
    if(step >= 3):            
        temperature_solver.Solve()
        fluid_solver.Solve()

    if(out == output_step): 
        PrintResults(model_part)  
        out = 0

    out = out + 1
    step = step + 1
      
gid_io.FinalizeResults()
          
        
