import fluid_ulf_var
##################################################################
##################################################################
#setting the domain size for the problem to be solved
      
def PrintSurfaceLevel(time, model_part, surfacelevelfilename):        
   y_max=-10.0
   for node in model_part.Nodes:
        if((node.X < 0.000001) and (node.GetSolutionStepValue(IS_STRUCTURE)==1) and (node.GetSolutionStepValue(IS_FLUID)==1)):
             y=node.Y
             if (y>y_max):
                  y_max=y
   outstring = str(time) + " " + str(y_max) + "\n"        
   surfacelevelfilename.write( outstring )


domain_size = 2

#importing applications
from KratosMultiphysics import *
from KratosMultiphysics.MeshingApplication import *
#from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ULFApplication import *
						                  


#defining a model part for the fluid and one for the structure
Qcomp_model_part = ModelPart("pFluidPart");

import runge_kutta_frac_step_solver
runge_kutta_frac_step_solver.AddVariables(Qcomp_model_part)

Qcomp_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
Qcomp_model_part.AddNodalSolutionStepVariable(VELOCITY)
Qcomp_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
Qcomp_model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
Qcomp_model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
Qcomp_model_part.AddNodalSolutionStepVariable(VISCOSITY)
Qcomp_model_part.AddNodalSolutionStepVariable(BODY_FORCE)
Qcomp_model_part.AddNodalSolutionStepVariable(IS_FLUID)
Qcomp_model_part.AddNodalSolutionStepVariable(VELOCITY_OLD)
Qcomp_model_part.AddNodalSolutionStepVariable(ADVPROJ)
Qcomp_model_part.AddNodalSolutionStepVariable(CONV_PROJ)
Qcomp_model_part.AddNodalSolutionStepVariable(PRESSUREAUX)
Qcomp_model_part.AddNodalSolutionStepVariable(VELOCITY_OLD_OLD)



input_file_name = fluid_ulf_var.problem_name
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_origin = ModelPartIO(input_file_name)
model_part_io_origin.ReadModelPart(Qcomp_model_part)

import math


Qcomp_model_part.SetBufferSize(3)

runge_kutta_frac_step_solver.AddDofs(Qcomp_model_part)


for node in Qcomp_model_part.Nodes:  
    node.SetSolutionStepValue(PRESSURE,0,0.0) 
    node.Free(VELOCITY_X) 
    node.Free(VELOCITY_Y) 
    node.Free(VELOCITY_Z) 
    node.Free(PRESSURE)
    if node.GetSolutionStepValue(IS_STRUCTURE)==1.0:
        node.Fix(VELOCITY_X)#
        node.Fix(VELOCITY_Y)
        node.Fix(VELOCITY_Z)
        node.SetSolutionStepValue(VELOCITY_X,0,0.0)
        node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
        node.SetSolutionStepValue(VELOCITY_Z,0,0.0)
    if node.GetSolutionStepValue(IS_FREE_SURFACE)==1.0:
        node.Fix(PRESSURE)
        node.SetSolutionStepValue(PRESSURE,0,0.0)


box_corner1 = Vector(3); 
box_corner1[0]=fluid_ulf_var.bounding_box_corner1_x; box_corner1[1]=fluid_ulf_var.bounding_box_corner1_y; box_corner1[2]=fluid_ulf_var.bounding_box_corner1_z;
box_corner2 = Vector(3); 
box_corner2[0]=fluid_ulf_var.bounding_box_corner2_x; box_corner2[1]=fluid_ulf_var.bounding_box_corner2_y; box_corner2[2]=fluid_ulf_var.bounding_box_corner2_z;


dynamic_tau=0.0
fluid_solver = runge_kutta_frac_step_solver.RungeKuttaFracStepSolver(Qcomp_model_part,domain_size, box_corner1, box_corner2)
fluid_solver.ReformDofAtEachIteration = False
pILUPrecond = ILU0Preconditioner()
fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-6, 5000,pILUPrecond)
fluid_solver.Initialize()


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
        
        res_name3 = str("free_surface_flow")
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
    
