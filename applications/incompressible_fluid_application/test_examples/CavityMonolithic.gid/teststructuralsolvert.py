##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2


from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *


#defining a model part
print "before creation of the model part"
model_part = ModelPart("FluidPart");
print "after creation of the model part"

##importing the solver files and adding the variables
import structural_solver_static
structural_solver_static.AddVariables(model_part)

#adding of Variables to Model Part should be here when the "very fix container will be ready"

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("cavityMonolithic",gid_mode,use_multifile,deformed_print_flag,write_conditions)
write_conditions = WriteConditionsFlag.WriteElementsOnly
##gid_io.ReadMesh(model_part.GetMesh())
gid_io.ReadModelPart(model_part)
##gid_io.WriteMesh((model_part).GetMesh(),domain_size,0.0,GiDPostMode.GiD_PostBinary);
gid_io.InitializeMesh( 0.0);
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

##add Degrees of Freedom to all of the nodes
structural_solver_static.AddDofs(model_part)



#creating a fluid solver object
fluid_solver = structural_solver_static.StaticStructuralSolver(model_part,domain_size)
#pILUPrecond = ILU0Preconditioner() 
#fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
##pDiagPrecond = DiagonalPreconditioner()
##fluid_solver.velocity_linear_solver =  BICGSTABSolver(1e-9, 5000,pDiagPrecond)
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pDiagPrecond)
##fluid_solver.pressure_linear_solver = SkylineLUFactorizationSolver();
fluid_solver.Initialize()

#settings to be changed
Re = 100
nsteps = 5
output_step = 1

Dt = 1000.0/Re
if(Dt > 0.1):
    Dt = 0.1
out = 0
for node in model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,1.0/Re)
    node.SetSolutionStepValue(BODY_FORCE_X,0,0)
    node.SetSolutionStepValue(BODY_FORCE_Y,0,0)

for node in model_part.Nodes:
    if(node.X < 0.001 and node.Y<0.001):
        node.Fix(PRESSURE)


zero = Vector(3);
zero[0] = 0.0;
zero[1] = 0.0;
zero[2] = 0.0;

gid_io.InitializeResults( 0.0, model_part.GetMesh() )

for step in range(1,nsteps):
    print "line49"

    time = Dt*step
    print time
    model_part.CloneTimeStep(time)

    print "qui"

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        fluid_solver.Solve()
        print "After solve"
        
    print "li"


    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(MESH_VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DENSITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VISCOSITY,model_part.Nodes,time,0)
        out = 0
    out = out + 1

node = model_part.Nodes[1]

gid_io.FinalizeResults();

print node

          
        

