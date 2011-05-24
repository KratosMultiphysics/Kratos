##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3

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
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosIncompressibleFluidApplication import *


#defining a model part
model_part = ModelPart("FluidPart");  

#providing the variable list to the model part
import fractional_step_solver
fractional_step_solver.AddVariables(model_part)
model_part.AddNodalSolutionStepVariable(REACTION)



#reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("cavity3D", gid_mode_flag, use_multifile, deformed_print_flag, write_conditions)
gid_io.ReadModelPart(model_part)

mesh_name = 0.0 #if we want the mesh to change at each time step then ****mesh_name = time****
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( model_part.GetMesh() )
gid_io.FinalizeMesh()

print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#adding dofs
fractional_step_solver.AddDofs(model_part)


#creating a fluid solver object
fluid_solver = fractional_step_solver.IncompressibleFluidSolver(model_part,domain_size)
#fluid_solver = flexible_incompressible_fluid_solver.IncompressibleFluidSolver(model_part,domain_size)
fluid_solver.laplacian_form =1;
fluid_solver.predictor_corrector = True;
fluid_solver.vel_toll = 1e-3
fluid_solver.time_order = 2
fluid_solver.ReformDofAtEachIteration = False
fluid_solver.echo_level = 0
fluid_solver.compute_reactions = True


##pILUPrecond = ILU0Preconditioner() 
##pDiagPrecond = DiagonalPreconditioner() 
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)
##fluid_solver.pressure_linear_solver = SkylineLUFactorizationSolver();

##for node in model_part.Nodes:
##    if(node.X < 0.001 and node.Y<0.001):
##        node.Fix(PRESSURE)

fluid_solver.Initialize()

outfile = open("reactions.out",'w')

def WriteReactions(time,nodes):
    outfile.write("time = " + str(time) + "\n" )
    for node in model_part.Nodes:
        vel = node.GetSolutionStepValue(REACTION)
        a = str(node.Id) + " " + str(vel[0]) + " " +  str(vel[1]) + " " + str(vel[2]) + "\n"
        outfile.write(a)
    outfile.write("\n")
    
#settings to be changed
Re = 100.0
nsteps = 100
output_step = 2

Dt = 1000.0/Re
if(Dt > 0.1):
    Dt = 0.1
out = 0
for node in model_part.Nodes:
    node.SetSolutionStepValue(DENSITY,0,1.0)
    node.SetSolutionStepValue(VISCOSITY,0,1.0/Re)

##zero = Vector(3);
##zero[0] = 0.0;
##zero[1] = 0.0;
##zero[2] = 0.0;
gid_io.InitializeResults(mesh_name , model_part.GetMesh())
for step in range(0,nsteps):
    print "line49"

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 4):
        fluid_solver.Solve()

        WriteReactions(time,model_part.Nodes)

##        for node in model_part.Nodes:
##            node.SetSolutionStepValue(PRESS_PROJ,0,zero);


    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
        out = 0
    out = out + 1
    
gid_io.FinalizeResults()

          
        

