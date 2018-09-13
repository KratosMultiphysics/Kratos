import fluid_only_var
#import fluid_ulf_var
##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

#importing applications
from KratosMultiphysics import *
from KratosMultiphysics.RadiationApplication import * 
from KratosMultiphysics.ExternalSolversApplication import *        #and now our application. note that we can import as many as we need to solve our specific
#defining a model part for the fluid and one for the structure


#setting the domain size for the problem to be solved
domain_size = 2  # 2D problem  

radiation_model_part = ModelPart("RadiationModelPart");

import Rad1


##RADIATION
SolverSettings_rad1 = Rad1.SolverSettings2
solver_constructor_rad1 = __import__(SolverSettings_rad1.solver_type)
solver_constructor_rad1.AddVariables(radiation_model_part, SolverSettings_rad1)

import nonlinear_convection_diffusionr_solver
nonlinear_convection_diffusionr_solver.AddVariables(radiation_model_part,SolverSettings_rad1) #the nodes are the same
print ("a44")


radiation_model_part.AddNodalSolutionStepVariable(INCIDENT_RADIATION_FUNCTION)
radiation_model_part.AddNodalSolutionStepVariable(RADIATIVE_INTENSITY_1)
radiation_model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
radiation_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
radiation_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)

#introducing input file name
input_file_name = fluid_only_var.problem_name
#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions 
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(radiation_model_part)

import math

radiation_model_part.SetBufferSize(3)

solver_constructor_rad1.AddDofs(radiation_model_part,SolverSettings_rad1)



radiation_model_part.Properties[0].SetValue(EMISSIVITY,1.0);
radiation_model_part.Properties[0].SetValue(AMBIENT_TEMPERATURE,298.0);
radiation_model_part.Properties[0].SetValue(CONVECTION_COEFFICIENT,0.0);
radiation_model_part.Properties[1].SetValue(EMISSIVITY,1.0);
radiation_model_part.Properties[1].SetValue(AMBIENT_TEMPERATURE,298.0);
radiation_model_part.Properties[1].SetValue(CONVECTION_COEFFICIENT,0.0);


#creating the solvers
radiation_solver1 = nonlinear_convection_diffusionr_solver.ConvectionDiffusionrSolver(radiation_model_part,domain_size,SolverSettings_rad1)
radiation_solver1.time_order = 1
radiation_solver1.ReformDofAtEachIteration = False
pPrecond = DiagonalPreconditioner()
pPrecond = DiagonalPreconditioner()
radiation_solver1.linear_solver = BICGSTABSolver(1e-6, 5000,pPrecond)
#radiation_solver1.linear_solver =  SkylineLUFactorizationSolver()
radiation_solver1.echo_level = 0
radiation_solver1.Initialize()

#mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( radiation_model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults(mesh_name,(radiation_model_part).GetMesh())

  
###
max_dt=0.002;
out = 0
time=0
step=0
output_step=0.000000010
next_output_time = output_step
final_time=1.0


#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<
#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<#>>>><<<<


for node in radiation_model_part.Nodes:
    node.SetSolutionStepValue(RADIATIVE_INTENSITY_1,0, 0.0)
    node.Free(RADIATIVE_INTENSITY_1)
    node.SetSolutionStepValue(INCIDENT_RADIATION_FUNCTION,0, 0.0)

for node in radiation_model_part.Nodes:
    node.SetSolutionStepValue(IS_BOUNDARY, 0, 0)
    node.SetSolutionStepValue(IS_INTERFACE, 0, 0)
    if(node.X<0.0001):
        node.SetSolutionStepValue(IS_BOUNDARY, 0, 1)
    if(node.X>0.999999):
        node.SetSolutionStepValue(IS_BOUNDARY, 0, 1)

for node in radiation_model_part.Nodes:
    node.SetSolutionStepValue(TEMPERATURE, 0, 500.0)
    if(node.X<0.0001):
        node.SetSolutionStepValue(TEMPERATURE, 0, 1000.0)
    if(node.X>0.999999):
        node.SetSolutionStepValue(TEMPERATURE, 0, 1000.0)

    
while(time < final_time):
   
              
    Dt = 0.0005; 

    time = time + Dt
    radiation_model_part.CloneTimeStep(time)

   
    radiation_solver1.Solve()

    for node in radiation_model_part.Nodes:
        TT = node.GetSolutionStepValue(RADIATIVE_INTENSITY_1)
        node.SetSolutionStepValue(INCIDENT_RADIATION_FUNCTION,0,TT)

    

    if(time > next_output_time):

        
        res_name4 = str("AIRE")
        gid_io.ChangeOutputName(res_name4)
        gid_io.InitializeMesh( time );
        gid_io.WriteNodeMesh((radiation_model_part).GetMesh());
        gid_io.WriteMesh((radiation_model_part).GetMesh());
        gid_io.FinalizeMesh();
        gid_io.InitializeResults(time, (radiation_model_part).GetMesh());
        gid_io.WriteNodalResults(TEMPERATURE,radiation_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_BOUNDARY,radiation_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(INCIDENT_RADIATION_FUNCTION,radiation_model_part.Nodes,time,0)
        gid_io.Flush()
        gid_io.FinalizeResults();



     
        
        next_output_time = next_output_time  + output_step;
	


        out = 0
    out = out + 1
    step = step + 1
    
