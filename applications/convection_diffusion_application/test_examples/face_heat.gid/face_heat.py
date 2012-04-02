##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

print "before importing kratos"
from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part

##########################################################
thermal_settings = ConvectionDiffusionSettings()
thermal_settings.SetDensityVariable(DENSITY)
thermal_settings.SetDiffusionVariable(CONDUCTIVITY)
thermal_settings.SetUnknownVariable(TEMPERATURE)
thermal_settings.SetVolumeSourceVariable(HEAT_FLUX)
thermal_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
thermal_settings.SetMeshVelocityVariable(MESH_VELOCITY)
thermal_settings.SetProjectionVariable(TEMP_CONV_PROJ);
##########################################################


import nonlinear_convection_diffusion_solver
nonlinear_convection_diffusion_solver.AddVariables(model_part,thermal_settings)

input_file_name = "face_heat"
#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#importing the solver files
nonlinear_convection_diffusion_solver.AddDofs(model_part,thermal_settings)

    
#creating a fluid solver object
solver = nonlinear_convection_diffusion_solver.ConvectionDiffusionSolver(model_part,domain_size,thermal_settings)
solver.time_order = 1
#pDiagPrecond = DiagonalPreconditioner()
#solver.linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)
solver.linear_solver = SkylineLUFactorizationSolver();
solver.echo_level = 0
solver.Initialize()

#assigning the fluid properties
conductivity = 0.25;
density = 900.0;
specific_heat = 2400.0;
temperature = 298.0;
for node in model_part.Nodes:
    node.SetSolutionStepValue(CONDUCTIVITY,0,conductivity);
    node.SetSolutionStepValue(DENSITY,0,density);
    node.SetSolutionStepValue(SPECIFIC_HEAT,0,specific_heat);
    node.SetSolutionStepValue(TEMPERATURE,0,temperature);

model_part.Properties[1][EMISSIVITY] = 1.0
model_part.Properties[1][AMBIENT_TEMPERATURE] = 25.0
model_part.Properties[1][CONVECTION_COEFFICIENT] = 8.0

   
#applying a temperature of 100
for node in model_part.Nodes:
    if(node.Y > 0.499):
         node.SetSolutionStepValue(FACE_HEAT_FLUX,0,1000.0);

 


#settings to be changed
nsteps = 10000
output_step = 50

Dt = 10.0

out = 1


for step in range(0,nsteps):
    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        solver.Solve()

    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(TEMP_CONV_PROJ,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(FACE_HEAT_FLUX,model_part.Nodes,time,0)
        out = 0
    out = out + 1

          
        

