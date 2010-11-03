##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs/' ##kratos_root/libs
kratos_applications_path = '../../../../applications/' ##kratos_root/applications
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_benchmarking_path)
#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_ULFApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosIncompressibleFluidApplication import *
from KratosULFApplication import *
import benchmarking


#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here when the "very fix container will be ready"
import runge_kutta_frac_step_solver
runge_kutta_frac_step_solver.AddVariables(model_part)

model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT)
model_part.AddNodalSolutionStepVariable(FORCE)

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("cilinderGLS",gid_mode,use_multifile,deformed_print_flag,write_conditions)
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


DragLift = Vector(2)

def BenchmarkCheck(time, model_part, section_nodes, center):
    
    DragLift=CalculateDragLift(time,section_nodes,center)
               
    benchmarking.Output(time, "Time")
    benchmarking.Output(DragLift[0], "drag", 0.0001)
    benchmarking.Output(DragLift[1], "lift", 0.0001)


##add Degrees of Freedom to all of the nodes
runge_kutta_frac_step_solver.AddDofs(model_part)

gravity = Vector(3);
gravity[0] = 0.0; gravity[1] = -1.0; gravity[2] = 0.0
zero = Vector(3);
zero[0] = 0.0; zero[1] = 0.0; zero[2] = 0.0
##for node in model_part.Nodes:    
##    if(node.Y > 0.99 and node.X > 0.001 and node.X < 0.999):
##        node.Free(VELOCITY_X);
##        node.Free(VELOCITY_Y);
##        #node.Free(VELOCITY_Z);
##    #node.Free(PRESSURE)
##    node.SetSolutionStepValue(BODY_FORCE,0,gravity);
##    node.SetSolutionStepValue(VELOCITY,0,zero);
        
#creating a fluid solver object
fluid_solver = runge_kutta_frac_step_solver.RungeKuttaFracStepSolver(model_part,domain_size)

##pILUPrecond = ILU0Preconditioner() 
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)

fluid_solver.Initialize()

#settings to be changed
Re = 100.0
nsteps = 300
output_step = 1

for node in model_part.Nodes :
    node.SetSolutionStepValue(VISCOSITY,0,1.0/Re);
    #node.SetSolutionStepValue(VISCOSITY,0,0.000017);
    node.SetSolutionStepValue(DENSITY,0,1.0);
    #node.SetSolutionStepValue(BODY_FORCE_Y,0,-10.00)

    #node.SetSolutionStepValue(PRESSURE,0,(1.0-node.Y)*10.0)
    #node.SetSolutionStepValue(VELOCITY_X,0,-1.0)
    #node.SetSolutionStepValue(VELOCITY_Y,0,-1.0)
    #node.SetSolutionStepValue(BODY_FORCE_Y,0,-50.0)
    #node.SetSolutionStepValue(BODY_FORCE_X,0,0.0)

#now we compute the delta time using CFL law
CFL_time_estimate_process=CFLProcess2D(model_part)
CFL=0.8;

#value to initialize only
Dt = 0.0;
dt_max=0.001;
out = 0
time=0.0

def SelectSectionNodes(section_nodes,model_part) :
    for node in model_part.Nodes:
        if (node.GetSolutionStepValue(IS_STRUCTURE)==1):
            section_nodes.append(node)
    print len(section_nodes)

#defining function to integrate forces on section contour
def CalculateCenterForces(time,section_nodes,center):
    
    fx = 0.0
    fy = 0.0
    mz = 0.0
    
    for node in section_nodes:
	fx -= node.GetSolutionStepValue(FORCE_X,0)
	fy -= node.GetSolutionStepValue(FORCE_Y,0)
	mz -=  + node.GetSolutionStepValue(FORCE_X,0) * (node.Y - center[1]) - node.GetSolutionStepValue(FORCE_Y,0) * (node.X - center[0])
    
    print "Drag =", fx
    print "Lift =", fy
    print "Moment =", mz
    
    return [fx,fy,mz]

def CalculateDragLift(time,section_nodes,center):
    
    fx = 0.0
    fy = 0.0
    mz = 0.0
    
    for node in section_nodes:
	fx -= node.GetSolutionStepValue(FORCE_X,0)
	fy -= node.GetSolutionStepValue(FORCE_Y,0)

    
    return [fx,fy] 
    
def PrintOutputFile(outputfile,time,Forces):
    
    outputfile.write(str(time)+" "+str(Forces[0])+" "+str(Forces[1])+" "+str(Forces[2])+"\n")
    
outstring2 = "Drag_Lift.txt"
outputfile = open(outstring2, 'w')

##################################################################
#################### Outputfile ###############################        
center = Vector(3)
center[0] = 0.0
center[1] = 0.0
center[2] = 0.0
#selecting section nodes
section_nodes = []
SelectSectionNodes(section_nodes,model_part);

print len(section_nodes);

gid_io.InitializeResults( 0.0, model_part.GetMesh() )

for step in range(0,nsteps):
    print "line49"
    Dt=(CFL_time_estimate_process).EstimateTime(CFL, dt_max)
    print "CFL gave this time step", Dt
           
    time = time + Dt
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        fluid_solver.Solve()
        Forces = CalculateCenterForces(time,section_nodes,center)
        #printing forces on section in the output file
        PrintOutputFile(outputfile,time,Forces)
        if (benchmarking.InBuildReferenceMode()):
            BenchmarkCheck(time, model_part, section_nodes, center)
        else:
            BenchmarkCheck(time, model_part, section_nodes, center)

    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(ACCELERATION,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(BODY_FORCE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_STRUCTURE,model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)

        out = 0
    out = out + 1

          
gid_io.FinalizeResults();
        

