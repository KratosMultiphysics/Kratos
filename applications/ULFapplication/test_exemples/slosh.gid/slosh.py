#domain size
domain_size = 2
#total simulation time
total_time = 10.0
#the max time step - it may be decreased in case it is necessary to avoid element inversion
max_delta_time = 0.01
#output time (every xxx seconds)
output_dt = 0.1
#safety factor for the delta time estimation
safety_factor = 0.5
#PATH of where the kratos library is installed
kratos_libs_path = '../../../../libs'
kratos_applications_path = '../../../../applications' 

project_name = 'slosh'
##################################################################
##################################################################
## ATTENTION: here the order is important

import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel


#importing applications
import applications_interface
applications_interface.Import_ULFApplication = True
applications_interface.Import_StructuralApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)


## from now on the order is not anymore crucial
##################################################################
##################################################################
print kernel
from KratosULFApplication import *
from KratosStructuralApplication import *

import math

def SelectSolidNodes(model_part,solid_nodes):
    print "in SelectSolidNodes"
    for node in model_part.Nodes:
        if(node.GetSolutionStepValue(IS_STRUCTURE) == 1):
            solid_nodes.append(node)   
    print solid_nodes

def FindNode(node_list,x,y,z):
    for node in node_list:
        a = (node.X - x)**2
        a = a + (node.Y - y)*(node.Y - y)
        a = a + (node.Z - z)*(node.Z - z)
        if (a < 0.00001):
            return node
    
def PrintForce(time,filename,node1,node2):
    print "in print force",time 
    outstring = str(time) + " "
    outstring += str(node1.GetSolutionStepValue(PRESSURE)) + " "
    outstring += str(node2.GetSolutionStepValue(PRESSURE)) + "\n"
    filename.write( outstring )            

        
def MoveSolidNodes(alpha_max_grad,T,Dt,time,solid_nodes,xc,yc):
    dist = Vector(3)
    out = Vector(3)
    disp = Vector(3)
    vel = Vector(3)
    center = Vector(3); center[0] = xc; center[1] = yc; center[2] = 0.0;
    
    #calculate rotation matrix
    alpha_rad = alpha_max_grad * 3.1415 / 180.0
    w = 2.0 * 3.1415 / T
    #incremental angle
    alpha = alpha_rad * ( math.sin(w * time) - math.sin( w * (time - Dt) ) );

    Rot = Matrix(3,3)
    Rot[0,0] = math.cos(alpha);  Rot[0,1] = -math.sin(alpha); Rot[0,2] = 0.0;
    Rot[1,0] = math.sin(alpha);  Rot[1,1] = math.cos(alpha);  Rot[1,2] = 0.0;
    Rot[2,0] = 0.0;              Rot[2,1] = 0.0;              Rot[2,2] = 0.0;    

        
    for i in range(0,len(solid_nodes)):
        dist[0] = solid_nodes[i].X 
        dist[1] = solid_nodes[i].Y 
        dist[2] = 0.0
        dist = dist - center
        
        #calcualte incremental displacement
        pos = Rot * dist + center

        disp[0] = pos[0] - solid_nodes[i].X;
        disp[1] = pos[1] - solid_nodes[i].Y;
        disp[2] = 0.0

        #calculate vel
        vel = disp / Dt
        solid_nodes[i].SetSolutionStepValue(VELOCITY,0,vel)

        #calcualte total displacement
        step_in_the_past = 1
        old_disp = solid_nodes[i].GetSolutionStepValue(DISPLACEMENT,step_in_the_past)
        disp = disp + old_disp
        solid_nodes[i].SetSolutionStepValue(DISPLACEMENT,0,disp)



#defining model parts
fluid_model_part = ModelPart("FluidPart");
structure_model_part = ModelPart("StructurePart");  
combined_model_part = ModelPart("CombinedPart");

#adding of Variables to Model Part should be here when the "very fix container will be ready"
#importing the solver files
import ulf_fsi
ulf_fsi.AddVariables(fluid_model_part)

#reading a model
#reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(project_name,gid_mode_flag,use_multifile,deformed_print_flag,write_conditions)
gid_io.ReadModelPart(fluid_model_part)
print fluid_model_part


#the buffer size should be set up here after the mesh is read for the first time
fluid_model_part.SetBufferSize(2)

ulf_fsi.AddDofs(fluid_model_part)

#setting the limits of the bounding box
box_corner1 = Vector(3); 
box_corner1[0]= -10.0; box_corner1[1]=-10.0; box_corner1[2]=-10.0;
box_corner2 = Vector(3); 
box_corner2[0]= 10.0; box_corner2[1]=10.0; box_corner2[2]=10.0;

#selecting the nodes for the motion of the sloshing
solid_nodes = []
SelectSolidNodes(fluid_model_part,solid_nodes)

#creating a fluid solver object
name = project_name
solver = ulf_fsi.ULF_FSISolver(fluid_model_part, structure_model_part, combined_model_part, box_corner1, box_corner2, domain_size)
solver.alpha_shape = 1.5;
solver.echo_level = 2;






Dt = 0.005
nsteps = 10000
#output_Dt = 0.05
output_Dt = output_dt
min_dt = 0.00001
max_dt = max_delta_time
safety_factor = 0.5 #you should put a safety factor ;-)!!!

next_output_time = output_Dt

#initializing the solver
solver.Initialize()


node_0_093 = FindNode(fluid_model_part.Nodes, 0.0, 0.093, 0.0)
node_0_299 = FindNode(fluid_model_part.Nodes, 0.0, 0.299, 0.0)

scalarout = open("pressure_history.out", 'w')
scalarout.write( "time (0.0,0.093) (0.0,0.299) \n")


time = 0.0
step = 0


for node in fluid_model_part.Nodes:
    press = (0.093 - node.Y) * 1000 * 9.81
    node.SetSolutionStepValue(PRESSURE,0,press);

while (time < total_time):   
    step = step + 1
    if(step<= 3):
	new_Dt = 0.000001
	new_Dt = safety_factor * new_Dt
    else:
    	new_Dt = solver.EstimateDeltaTime(max_dt,domain_size)
    	new_Dt = safety_factor * new_Dt
    #print "forever dt", new_Dt
    #0.5 - SAFETY FACTOR
##    if (new_Dt==Dt):
    time = time + new_Dt
##    else:
##        time = time + 0.1*new_Dt
    
####    time = Dt*step
##    new_Dt = Dt
    fluid_model_part.CloneTimeStep(time)
    structure_model_part.CloneTimeStep(time)
    combined_model_part.CloneTimeStep(time)

    #moving the gate
    alpha_max_grad = 4.0;
    T = 1.91;
    xc = 0.45; yc=0.184;
    MoveSolidNodes(alpha_max_grad,T,new_Dt,time,solid_nodes,xc,yc)


    #solving the fluid problem
    if(step > 3):
        
        print "before the solution"

        solver.Solve()
        
        print "after completing the solution"

        if(time > next_output_time):
    
            file_name = project_name
            file_name = file_name + str(time)
            gid_io.InitializeMesh( time );
            gid_io.WriteNodeMesh((combined_model_part).GetMesh());
            gid_io.WriteMesh((combined_model_part).GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(time, (combined_model_part).GetMesh());
            
            gid_io.WriteNodalResults(DISPLACEMENT, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(NODAL_H, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_FLUID, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_BOUNDARY, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_FREE_SURFACE, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_STRUCTURE, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(VELOCITY, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(PRESSURE, (combined_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(BODY_FORCE, (combined_model_part).Nodes, time, 0);
            gid_io.FinalizeResults();

            next_output_time = next_output_time  + output_Dt;
          
        


