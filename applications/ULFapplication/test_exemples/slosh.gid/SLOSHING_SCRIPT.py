###################################################################################FUNCTIONS FOR MOVING THE CONTAINER####################################################
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


#####################################################################################################################################################################

import fluid_ulf_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = fluid_ulf_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = fluid_ulf_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = fluid_ulf_var.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_ULFApplication = True
applications_interface.Import_MeshingApplication = True
applications_interface.Import_PFEMApplication = True
applications_interface.Import_StructuralApplication = True
##applications_interface.Import_ExternalSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosULFApplication import *
from KratosMeshingApplication import *
from KratosStructuralApplication import *

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");
structure_model_part = ModelPart("StructurePart");  
combined_model_part = ModelPart("CombinedPart");
fluid_only_model_part = ModelPart("FluidOnlyPart");



fluid_model_part.AddNodalSolutionStepVariable(IS_VISITED);
fluid_model_part.AddNodalSolutionStepVariable(DISTANCE);

SolverType=fluid_ulf_var.SolverType
#if (SolverType=="Incompressible_Modified_FracStep" or SolverType=="FracStep"):
    #fluid_only_model_part = ModelPart("FluidOnlyPart");

#############################################
##importing the solvers needed
if(SolverType == "Incompressible_Modified_FracStep"):
    import ulf_frac
    ulf_frac.AddVariables(fluid_model_part)
    print "You are using a modified ULF FRAC solver------------------------------------------------------------------------------------------------------------------"
elif(SolverType == "FracStep"):
    import ulf_frac
    ulf_frac.AddVariables(fluid_model_part)       
elif(SolverType == "Quasi_Inc_Constant_Pressure"):
    import ulf_fsi
    ulf_fsi.AddVariables(fluid_model_part)
elif(SolverType == "Quasi_Inc_Linear_Pressure"):
    import ulf_fsi_inc
    ulf_fsi_inc.AddVariables(fluid_model_part)
else:
    raise "solver type not supported: options are fractional_step - \
	modified_frac_steop - quasi_inc_constant_pres - \
	quasi_inc_lin_pres"

#introducing input file name
input_file_name = fluid_ulf_var.problem_name

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_origin = ModelPartIO(input_file_name)
model_part_io_origin.ReadModelPart(fluid_model_part)



#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)
combined_model_part.SetBufferSize(3)
compute_reactions=1

##adding dofs
if(SolverType == "Incompressible_Modified_FracStep"):
    ulf_frac.AddDofs(fluid_model_part, compute_reactions)
elif(SolverType == "FracStep"):
    ulf_frac.AddDofs(fluid_model_part, compute_reactions)
elif(SolverType == "Quasi_Inc_Constant_Pressure"):
    ulf_fsi.AddDofs(fluid_model_part, compute_reactions)
elif(SolverType == "Quasi_Inc_Linear_Pressure"):
    ulf_fsi_inc.AddDofs(fluid_model_part, compute_reactions)


#if(SolverType == "Quasi_Inc_Constant_Pressure" or SolverType == "Quasi_Inc_Linear_Pressure"):
#      for node in fluid_model_part.Nodes:
#	  node.Free(PRESSURE)

#setting the limits of the bounding box
box_corner1 = Vector(3); 
box_corner1[0]=fluid_ulf_var.bounding_box_corner1_x; box_corner1[1]=fluid_ulf_var.bounding_box_corner1_y; box_corner1[2]=fluid_ulf_var.bounding_box_corner1_z;
box_corner2 = Vector(3); 
box_corner2[0]=fluid_ulf_var.bounding_box_corner2_x; box_corner2[1]=fluid_ulf_var.bounding_box_corner2_y; box_corner2[2]=fluid_ulf_var.bounding_box_corner2_z;

#here we write the convergence data..,
outstring2 = "convergence_info.txt"
outputfile1 = open(outstring2, 'w')

#selecting the nodes for the motion of the sloshing
solid_nodes = []
SelectSolidNodes(fluid_model_part,solid_nodes)

add_nodes=fluid_ulf_var.adaptive_refinement
bulk_modulus=fluid_ulf_var.bulk_modulus
density=fluid_ulf_var.density
#creating the solvers
#fluid solver

FSI=0

if(SolverType == "Incompressible_Modified_FracStep"):        
    solver = ulf_frac.ULF_FSISolver(outputfile1, fluid_only_model_part, fluid_model_part, structure_model_part, combined_model_part, FSI, compute_reactions, box_corner1, box_corner2, domain_size, add_nodes, bulk_modulus, density)
    solver.alpha_shape = 1.5#fluid_ulf_var.alpha_shape;
    solver.echo_level = 2;
    
    for node in fluid_model_part.Nodes:
	node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)
	node.SetSolutionStepValue(DENSITY,0, density)
	node.SetSolutionStepValue(VISCOSITY,0, 0.0001)
	node.SetSolutionStepValue(BODY_FORCE_Y,0, -10.000)
    solver.Initialize()
    
if(SolverType == "FracStep"):    
    solver = ulf_frac.ULF_FSISolver(outputfile1, fluid_only_model_part, fluid_model_part, structure_model_part, combined_model_part, compute_reactions, box_corner1, box_corner2, domain_size, add_nodes, bulk_modulus, density)
    solver.alpha_shape = fluid_ulf_var.alpha_shape;
    solver.echo_level = 2;
    for node in fluid_model_part.Nodes:
	node.SetSolutionStepValue(BULK_MODULUS,0, 0.0)   
	node.SetSolutionStepValue(DENSITY,0, density)   
    solver.Initialize()
    
elif(SolverType == "Quasi_Inc_Constant_Pressure"):
    solver = ulf_fsi.ULF_FSISolver(fluid_model_part, structure_model_part, combined_model_part, compute_reactions, box_corner1, box_corner2, domain_size, add_nodes)
    solver.alpha_shape = fluid_ulf_var.alpha_shape;
    solver.echo_level = 2;
   
    
    for node in fluid_model_part.Nodes:
	node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)   
	node.SetSolutionStepValue(DENSITY,0, density)   
    solver.Initialize()
    
elif(SolverType == "Quasi_Inc_Linear_Pressure"): 
    solver = ulf_fsi_inc.ULF_FSISolver(outputfile1, fluid_model_part, structure_model_part, combined_model_part, compute_reactions,  box_corner1, box_corner2, domain_size, add_nodes, bulk_modulus, density)
    solver.alpha_shape = fluid_ulf_var.alpha_shape;
    solver.echo_level = 2;
        
    for node in fluid_model_part.Nodes:
	node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)   
	node.SetSolutionStepValue(DENSITY,0, density)   
    solver.Initialize()

##check to ensure that no node has zero density or pressure
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise 'node with zero VISCOSITY found'



print "fluid solver created"

#settings to be changed
Dt = fluid_ulf_var.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01
final_time = fluid_ulf_var.max_time
output_step = fluid_ulf_var.output_step
safety_factor = 0.5 #you should put a safety factor ;-)!!!

next_output_time = output_step

time = 0.0
step = 0

inlet_vel = Vector(3)
inlet_vel[0]=0.0
inlet_vel[1]=0.0
inlet_vel[2]=0.0


def SelectVisited(nodes):
     for node in nodes:
         if(node.GetSolutionStepValue(IS_FREE_SURFACE)==1.0):
             node.SetValue(IS_VISITED,1.0)
             node.SetSolutionStepValue(DISTANCE,0,0.0)
             #print "AAAAAAAAAAAAA"
         else:
             node.SetValue(IS_VISITED,0.0)

#set_h_map_process = SetHMapProcess(fluid_model_part);

SelectVisited(fluid_model_part.Nodes)

distance_utils=BodyDistanceCalculationUtils()

dummy=LagrangianInletProcess(fluid_model_part, 0.0, inlet_vel)

while (time < final_time):
    step = step+1   
    
    
    print time
    if(step <= 3):
        new_Dt = 0.00000001;
        time = time + new_Dt*safety_factor

    #solving the fluid problem
    if(step > 3):
        new_Dt = solver.EstimateDeltaTime(Dt, domain_size)
        time = time + new_Dt*safety_factor
        
        
        combined_model_part.CloneTimeStep(time)
        
	alpha_max_grad = 4.0;
	T = 1.91;
	xc = 0.45; yc=0.184;
	MoveSolidNodes(alpha_max_grad,T,new_Dt,time,solid_nodes,xc,yc)

        solver.Solve(dummy)        
        print "after completing the solution"
##        SelectVisited(fluid_model_part.Nodes)
##        distance_utils.CalculateDistances2D(fluid_model_part.Elements, DISTANCE, True)
##
##        for node in fluid_model_part.Nodes:
##            if (node.GetSolutionStepValue(DISTANCE)<=0.09 and node.GetSolutionStepValue(DISTANCE)>=0.0001):
##                h=0.05*node.GetSolutionStepValue(DISTANCE)+0.003
##                node.SetSolutionStepValue(NODAL_H, 0, h)
##            else:
##                node.SetSolutionStepValue(NODAL_H, 0, 0.01)
                  

##        if (time>0.05):
##            #here we want to store the min_H, max_H, min_dist, max_dist
##            min_H=1000000.0
##            max_H=0.0
##            min_dist=0.0
##            max_dist=0.0
##            for node in fluid_model_part.Nodes:
##                H=node.GetSolutionStepValue(NODAL_H)
##                if (min_H>H):
##                    min_H=H
##                if (max_H<H):
##                    max_H=H
##
##            #min_H=1.0*min_H
##            #we set min_Dist to the min_H
##            #I want to prescribe the mesh size: min_H=0.01
##            min_H=0.004
##            min_Dist=0.5*min_H
##            #1.0*min_H
##
##            set_h_map_process.CalculateOptimalH(min_H, max_H)

        if(time > next_output_time):
    
            file_name = input_file_name
            file_name = file_name + str(time)

            gid_io.InitializeMesh( time );
            gid_io.WriteNodeMesh((combined_model_part).GetMesh());
            gid_io.WriteMesh((combined_model_part).GetMesh());
            gid_io.FinalizeMesh();

            gid_io.InitializeResults(time, (combined_model_part).GetMesh());

            gid_io.WriteNodalResults(DISTANCE, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(NODAL_H, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_VISITED, combined_model_part.Nodes, time, 0);           
            gid_io.WriteNodalResults(IS_FREE_SURFACE, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_INTERFACE, combined_model_part.Nodes, time, 0);            
            gid_io.WriteNodalResults(VELOCITY, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(PRESSURE, (combined_model_part).Nodes, time, 0);
            
            

            gid_io.Flush()
            #gid_io.CloseResultFile();
            gid_io.FinalizeResults()

            next_output_time = next_output_time  + output_step;
 
