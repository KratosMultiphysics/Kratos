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

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_ThermoMechanicalApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)


#defining a model part
model_part = ModelPart("FluidPart");  


## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosThermoMechanicalApplication import *

##########################################################
thermal_settings = ConvectionDiffusionSettings()
thermal_settings.SetDensityVariable(DENSITY)
thermal_settings.SetDiffusionVariable(CONDUCTIVITY)
thermal_settings.SetUnknownVariable(TEMPERATURE)
thermal_settings.SetVolumeSourceVariable(HEAT_FLUX)
thermal_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
thermal_settings.SetMeshVelocityVariable(MESH_VELOCITY)
##########################################################

#importing the solver files
import thermo_monolithic_solver_eulerian
thermo_monolithic_solver_eulerian.AddVariables(model_part,thermal_settings)




#adding of Variables to Model Part should be here when the "very fix container will be ready"
#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("square_split",gid_mode,multifile,deformed_mesh_flag, write_conditions)
gid_io.ReadModelPart(model_part)

#mesh_name = 0.0
#gid_io.InitializeMesh( mesh_name );
#gid_io.WriteMesh((model_part).GetMesh());
#gid_io.FinalizeMesh()
#print model_part


#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#importing the solver files
thermo_monolithic_solver_eulerian.AddDofs(model_part)
print "111111111111111111111111111111111111111111"

#creating a fluid solver object
solver = thermo_monolithic_solver_eulerian.MonolithicSolver(model_part,domain_size,thermal_settings)
a = Matrix(1,2)
a[0,0] = 1
a[0,1] = 2
solver.contact_matix =  a
solver.Namx = len(model_part.Properties)
#pDiagPrecond = DiagonalPreconditioner()
#solver.linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
solver.Initialize()

#assigning the fluid properties
conductivity = 0.10;
density = 1.0;
specific_heat = 1.0;
for node in model_part.Nodes:
    node.SetSolutionStepValue(CONDUCTIVITY,0,conductivity);
    node.SetSolutionStepValue(DENSITY,0,density);
    node.SetSolutionStepValue(SPECIFIC_HEAT,0,specific_heat);

#assigning a rotational velocity field
for node in model_part.Nodes:
    node.SetSolutionStepValue(TEMPERATURE,0,0.0);
    if(node.Y > 0.499 ):
         node.SetSolutionStepValue(TEMPERATURE,0,5.0);
         node.Fix(TEMPERATURE)

#for node in model_part.Nodes:
    #node.SetSolutionStepValue(VELOCITY_Y,0,-0.10);
    #node.Fix(VELOCITY_Y)  

#settings to be changed
nsteps = 6000
output_step = 100

Dt = 0.0020;

out = 10

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())



for step in range(0,nsteps):
    print "line49"
    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]    #solving the fluid problem
    if(step > 3):
     solver.Solve()

    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(NODAL_AREA,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(SPECIFIC_HEAT,model_part.Nodes,time,0) 
        gid_io.WriteNodalResults(DENSITY,model_part.Nodes,time,0) 
        gid_io.WriteNodalResults(CONDUCTIVITY,model_part.Nodes,time,0)      
        gid_io.WriteNodalResults(MESH_VELOCITY,model_part.Nodes,time,0)    
        gid_io.WriteNodalResults(NORMAL,model_part.Nodes,time,0)           
        
        out = 0
    out = out + 1
gid_io.FinalizeResults()

          
        

