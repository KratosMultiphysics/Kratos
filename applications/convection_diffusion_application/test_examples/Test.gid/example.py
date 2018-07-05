import fluid_ulf_var

domain_size = 2

#importing applications
import sys
from KratosMultiphysics import *
#from KratosMultiphysics.MeshingApplication import *
#from KratosMultiphysics.ExternalSolversApplication import *
#from KratosMultiphysics.ULFApplication import *
#sys.path.append(ProjectParameters.kratos_path)
#from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
#from KratosMultiphysics.ExternalSolversApplication import *
#from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
#defining a model part for the fluid and one for the structure
model_part = ModelPart("pFluidPart");

#import runge_kutta_frac_step_solver
#runge_kutta_frac_step_solver.AddVariables(Qcomp_model_part)

import coupled_fluid_thermal_solver

#coupled_fluid_thermal_solver.AddVariables(Qcomp_model_part)

#Qcomp_model_part.AddNodalSolutionStepVariable(VELOCITY)
#Qcomp_model_part.AddNodalSolutionStepVariable(PRESSURE)



input_file_name = "Test"
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_origin = ModelPartIO(input_file_name)
model_part_io_origin.ReadModelPart(model_part)


#parameters = "ProjectParameters.json"
with open("ProjectParameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

#SSSSSSS

#Qcomp_model_part.SetBufferSize(3)

#runge_kutta_frac_step_solver.AddDofs(Qcomp_model_part)

#with open("ProjectParameters.json",'r') as parameter_file:
#    parameters = KratosMultiphysics.Parameters(parameter_file.read())
#ddddddddddd
domain_size=3
#dynamic_tau=0.0

fluid_solver = coupled_fluid_thermal_solver.CoupledFluidThermalSolver(model_part,parameters)
ssssssssssssssss
fluid_solver.ReformDofAtEachIteration = False
pILUPrecond = ILU0Preconditioner()
fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-6, 5000,pILUPrecond)
fluid_solver.Initialize()

ssssssssss
#mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( Qcomp_model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults(mesh_name,(Qcomp_model_part).GetMesh())

   
###TIME
out = 0
time=0
step=0
output_step=0.000001
next_output_time = output_step
final_time=1.0

###
domain_size=2

ulf_time_step_dec_process = UlfTimeStepDecProcess(Qcomp_model_part);


for node in Qcomp_model_part.Nodes :
    node.SetSolutionStepValue(BODY_FORCE_Y,0,-10.00)
    node.SetSolutionStepValue(VISCOSITY,0,0.00001);
    node.SetSolutionStepValue(DENSITY,0,1000.0);

outstring1 = "sensor_height.txt"
outputfile1 = open(outstring1, 'w')
out_file = outputfile1
out_file.write("WAVE HEIGHT AT SENSOR 1" + "\n")

while(time < final_time):
   

    Dt=0.005

    time=Qcomp_model_part.ProcessInfo[TIME]

    time = time + Dt

    Qcomp_model_part.CloneTimeStep(time)
    
    
    fluid_solver.Solve()
    PrintSurfaceLevel(time, Qcomp_model_part, outputfile1)

               
    
    if(time > next_output_time):
        
        res_name3 = str("POLIMERO")
        gid_io.ChangeOutputName(res_name3)
        gid_io.InitializeMesh( time );
        gid_io.WriteNodeMesh((Qcomp_model_part).GetMesh());
        gid_io.WriteMesh((Qcomp_model_part).GetMesh());
        gid_io.FinalizeMesh();
        gid_io.InitializeResults(time, (Qcomp_model_part).GetMesh());
        gid_io.WriteNodalResults(VELOCITY,Qcomp_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,Qcomp_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_FREE_SURFACE,Qcomp_model_part.Nodes,time,0)
        gid_io.Flush()
        gid_io.FinalizeResults();
        
        next_output_time = next_output_time  + output_step;
	

        out = 0
    out = out + 1
    step = step + 1
    
