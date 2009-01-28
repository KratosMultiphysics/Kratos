#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *
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
        if (a < 0.0000001):
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
        
            
            

#path for applications s
import sys
sys.path.append('/home/rrossi/kratos_cvs/applications/incompressible_fluid_application/python_scripts')

kernel = Kernel()

#############################################
#defining the domain size
domain_size = 2

#############################################
#############################################
#attention!! order is crucial here!!!
#importing the applications library and adding them to the kernel
incompressible_fluid_application = KratosIncompressibleFluidApplication()
kernel.AddApplication(incompressible_fluid_application)

#dynamic renumbering of variables to ensure consistency 
kernel.Initialize()
kernel.InitializeApplication(incompressible_fluid_application);
#############################################
#############################################

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here when the "very fix container will be ready"

#reading a model
gid_io = GidIO("sloshing_experimental",GiDPostMode.GiD_PostBinary)
gid_io.ReadMesh(model_part.GetMesh())
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
import pfem_solver
pfem_solver.AddVariables(model_part)

#setting the limits of the bounding box
box_corner1 = Vector(3); box_corner1[0]=-1.5; box_corner1[1]=-1.5; box_corner1[2]=-0.5;
box_corner2 = Vector(3); box_corner2[0]=1.5; box_corner2[1]=1.5;  box_corner2[2]=0.1;

#selecting the nodes for the motion of the gate
solid_nodes = []
SelectSolidNodes(model_part,solid_nodes)

#creating a fluid solver object
name = str("sloshing_experimental")

##solver = pfem_solver.PFEMSolver(model_part,name,box_corner1,box_corner2,domain_size)
##pDiagPrecond = DiagonalPreconditioner()
##solver.pressure_linear_solver =  CGSolver(1e-3, 5000,pDiagPrecond)
#solver.predictor_corrector = True
Mesher = TriGenModeler()

Dt = 0.005
nsteps = 6
output_Dt = 0.05
min_dt = 0.2*Dt
max_dt = 2.0*Dt
tmax = 0.35


node_0_093 = FindNode(model_part.Nodes, 0.0, 0.093, 0.0)
node_0_299 = FindNode(model_part.Nodes, 0.0, 0.299, 0.0)

scalarout = open("pressure_history.out", 'w')
scalarout.write( "time (0.0,0.093) (0.0,0.299) \n")

time = 0.0
step = 0
#for step in range(0,nsteps):
while (time < tmax):
    step = step + 1
    new_Dt = Dt

    time = time + new_Dt
    print "time = ", time, " new_Dt= ",new_Dt," step = ", step 

    #time = Dt*step
    model_part.CloneTimeStep(time)

    #moving the gate
    alpha_max_grad = 4.0;
    T = 1.91;
    xc = 0.45; yc=0.184;
    MoveSolidNodes(alpha_max_grad,T,new_Dt,time,solid_nodes,xc,yc)
    
    print time
    (Mesher).ReGenerateMesh(model_part,1.5)
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    





          
        

