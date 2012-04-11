#setting the domain size for the problem to be solved
domain_size = 2

#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *

#importing system libraries as needed
import math


def SelectNodes(moving_nodes,model_part):
    for node in model_part.Nodes:
        if(node.Y > 0.49 and node.Y < 0.51):
            if(node.X > 0.0001 and node.X < 0.99999 ):
                moving_nodes.append(node)
                node.Fix(DISPLACEMENT_X);
                node.Fix(DISPLACEMENT_Y);
                node.Fix(DISPLACEMENT_Z);
                
def ApplyDisplacementConditions(model_part):
    for node in model_part.Nodes:
        if(node.Y > 0.9999 or node.Y < 0.0001 or node.X < 0.0001 or node.X > 0.999):
            node.Fix(DISPLACEMENT_X);
            node.Fix(DISPLACEMENT_Y);
            node.Fix(DISPLACEMENT_Z);
                
def MoveNodes(moving_nodes,time):
    disp = Vector(3)
    disp[0] = 0.0; disp[2] = 0.0
    disp[1] = 0.05*math.sin(time*1.0)
    print disp
    for node in moving_nodes:
        node.SetSolutionStepValue(DISPLACEMENT,0,disp);
        
#defining a model part
model_part = ModelPart("FluidPart");  

#importing the solver files
import incompressible_fluid_solver
print "qui"
import mesh_solver
print "li"

#adding of Variables to Model Part should be here when the "very fix container will be ready"
incompressible_fluid_solver.AddVariables(model_part)
mesh_solver.AddVariables(model_part)

#reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("cavity2d", gid_mode_flag, use_multifile, deformed_print_flag, write_conditions)

gid_io.ReadModelPart(model_part)

mesh_name = 0.0 #if we want the mesh to change at each time step then ****mesh_name = time****
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( model_part.GetMesh() )
gid_io.FinalizeMesh()
        
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#adding Dofs
incompressible_fluid_solver.AddDofs(model_part)
mesh_solver.AddDofs(model_part)


#creating a fluid solver object
fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(model_part,domain_size)
fluid_solver.laplacian_form =1;
fluid_solver.vel_toll = 1e-3
fluid_solver.time_order = 2

##pILUPrecond = ILU0Preconditioner() 
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
fluid_solver.Initialize()

#creating a mesh solver
print "qui"
import mesh_solver
reform_dofs_at_each_step = False
mesh_sol = mesh_solver.MeshSolver(model_part,domain_size,reform_dofs_at_each_step)
mesh_sol.Initialize()

#settings to be changed
Re = 10000.0
nsteps = 3000
output_step = 10

Dt = 1000.0/Re
if(Dt > 0.1):
    Dt = 0.1
out = 0
for node in model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,1.0/Re)


#selecting nodes to be moved
moving_nodes = []
SelectNodes(moving_nodes,model_part);
ApplyDisplacementConditions(model_part);
print len(moving_nodes);



gid_io.InitializeResults(mesh_name , model_part.GetMesh()) #gid_io.InitializeResults(time , model_part.GetMesh()) if gid_io.InitializeMesh( time )

#solve
for step in range(0,nsteps):
    print "line49"

    

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    
    
    #solving the fluid problem
    if(step > 3):
        MoveNodes(moving_nodes,time)
        mesh_sol.Solve();
        
        fluid_solver.Solve()


    #print the results
    if(out == output_step):
        #results writing options
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(MESH_VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
        gid_io.Flush()
        out = 0
    out = out + 1

gid_io.FinalizeResults()
         
        

