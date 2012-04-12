##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.PFEMApplication import *

###########################################################
thermal_settings = ConvectionDiffusionSettings()
thermal_settings.SetDensityVariable(DENSITY)
thermal_settings.SetDiffusionVariable(CONDUCTIVITY)
thermal_settings.SetUnknownVariable(TEMPERATURE)
thermal_settings.SetVolumeSourceVariable(HEAT_FLUX)
thermal_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
thermal_settings.SetMeshVelocityVariable(MESH_VELOCITY)
##########################################################
from KratosPFEMApplication import *

#defining a model part
model_part = ModelPart("FluidPart");  
temperature_model_part = ModelPart("TemperaturePart");

#adding of Variables to Model Part should be here when the "very fix container will be ready"
import runge_kutta_frac_step_comp_solver
runge_kutta_frac_step_comp_solver.AddVariables(model_part)

model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT)
model_part.AddNodalSolutionStepVariable(FORCE)
model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
model_part.AddNodalSolutionStepVariable(TEMPERATURE)

import nonlinear_convection_diffusion_solver
nonlinear_convection_diffusion_solver.AddVariables(model_part, thermal_settings) #the nodes are the same 


#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("mixing",gid_mode,use_multifile,deformed_print_flag,write_conditions)
write_conditions = WriteConditionsFlag.WriteElementsOnly
#gid_io.ReadMesh(model_part.GetMesh())
gid_io.ReadModelPart(model_part)
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh();
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)


#for node in model_part.Nodes:
#    if (node.X<0.0001):
    #if (node.X>0.2 and node.X<0.5 and node.Y<0.00001):
#        node.SetSolutionStepValue(FLAG_VARIABLE,0,3.0)

##add Degrees of Freedom to all of the nodes
runge_kutta_frac_step_comp_solver.AddDofs(model_part)
nonlinear_convection_diffusion_solver.AddDofs(model_part, thermal_settings)

NISTTools = NistUtils()

gravity = Vector(3);
gravity[0] = 0.0; gravity[1] = -1.0; gravity[2] = 0.0
zero = Vector(3);
zero[0] = 0.0; zero[1] = 0.0; zero[2] = 0.0
for node in model_part.Nodes:
    node.Free(PRESSURE)
    
for node in model_part.Nodes:
    if (node.Y<0.0001 and node.X<0.0001):
        node.SetSolutionStepValue(PRESSURE,0,0.0)
        node.Fix(PRESSURE)
    
#    if(node.Y > 0.99):
#        node.SetSolutionStepValue(PRESSURE, 0, 0.0)
#        node.Fix(PRESSURE)
#        node.Free(VELOCITY_X);
#        node.Free(VELOCITY_Y);
##        #node.Free(VELOCITY_Z);
    
##    node.SetSolutionStepValue(BODY_FORCE,0,gravity);
##    node.SetSolutionStepValue(VELOCITY,0,zero);
        
#creating a fluid solver object
fluid_solver = runge_kutta_frac_step_comp_solver.RungeKuttaFracStepCompSolver(model_part,domain_size)

##pILUPrecond = ILU0Preconditioner() 
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)

fluid_solver.Initialize()

#convection diffusion solver
temperature_solver = nonlinear_convection_diffusion_solver.ConvectionDiffusionSolver(temperature_model_part,domain_size, thermal_settings)
temperature_solver.time_order = 1
temperature_solver.ReformDofAtEachIteration = False
temperature_solver.echo_level = 0
temperature_solver.Initialize()

utilities = VariableUtils()

#settings to be changed
Re = 10.0
nsteps = 10000
output_step = 20

for node in model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,0.001)
    node.SetSolutionStepValue(DENSITY,0,1.0000)
    node.SetSolutionStepValue(TEMPERATURE,0,293.16)
    #node.SetSolutionStepValue(PRESSURE,0,(1.0-node.Y)*1.0)
    #node.SetSolutionStepValue(VELOCITY_X,0,-1.0)
    node.SetSolutionStepValue(BODY_FORCE_Y,0,-10.0)
    #prescribe hydrostatic pressure
    node.SetSolutionStepValue(PRESSURE,0,1.0*(1.0-node.Y) )


#for node in model_part.Nodes:
#    node.Free(PRESSURE)
#    if(node.Y > 0.99):
#        node.Fix(PRESSURE)

conductivity = 1.0;

specific_heat = 1000.0;

for node in model_part.Nodes:
    node.SetSolutionStepValue(CONDUCTIVITY,0,conductivity);    
    node.SetSolutionStepValue(SPECIFIC_HEAT,0,specific_heat);

model_part.Properties[1][EMISSIVITY] = 0.0
model_part.Properties[1][AMBIENT_TEMPERATURE] = 293.16
model_part.Properties[1][CONVECTION_COEFFICIENT] = 0.0

temperature_neigh_finder = FindNodalNeighboursProcess(temperature_model_part,9,18)

##applying temperature boundary conditions
for node in model_part.Nodes:
    if (node.X<0.0001):
    #if (node.X>0.2 and node.X<0.5 and node.Y<0.00001):
        node.SetSolutionStepValue(TEMPERATURE,0,298.16)
        node.Fix(TEMPERATURE)
    if (node.X>0.9999):
        node.SetSolutionStepValue(TEMPERATURE,0,288.16)
        node.Fix(TEMPERATURE)

#now we compute the delta time using CFL law
CFL_time_estimate_process=CFLProcess2D(model_part)
CFL=0.8;

dt_max = 0.01
Dt=0.1
out = 0
time=0.0

gid_io.InitializeResults( 0.0, model_part.GetMesh() )

for step in range(0,nsteps):
    print "line49"
    Dt=(CFL_time_estimate_process).EstimateTime(CFL, dt_max)
    print "CFL gave this time step", Dt
    time = time + Dt

    model_part.CloneTimeStep(time)
    temperature_model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        NISTTools.GenerateModelPart(model_part,temperature_model_part,domain_size);
        temperature_solver.Solve()
##        for node in model_part.Nodes:
##            grav=-10.0+10.0*(0.01*( node.GetSolutionStepValue(TEMPERATURE)-293.16))
##            node.SetSolutionStepValue(BODY_FORCE_Y,0,grav); 
        fluid_solver.Solve()
        #temperature_solver.Solve()
        (temperature_model_part.Nodes).clear()
        (temperature_model_part.Elements).clear()
        (temperature_model_part.Conditions).clear()
        (temperature_model_part.Properties).clear()

    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DENSITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(ACCELERATION,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(BODY_FORCE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)

        out = 0
    out = out + 1

          
gid_io.FinalizeResults();
        

