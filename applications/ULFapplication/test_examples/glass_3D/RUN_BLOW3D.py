#import problem_settings
import problem_settings

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = problem_settings.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
import sys
sys.path.append(problem_settings.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.StructuralApplication import *

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");
structure_model_part = ModelPart("StructurePart");  
combined_model_part = ModelPart("CombinedPart");
wall_model_part = ModelPart("WallPart");


SolverType=problem_settings.SolverType
if (SolverType=="Incompressible_Modified_FracStep" or SolverType=="FracStep"):
    fluid_only_model_part = ModelPart("FluidOnlyPart");

#wall_model_part = ModelPart("WallPart");
#############################################
##importing the solvers needed
if(SolverType == "Incompressible_Modified_FracStep"):
    import ulf_PGLASS
    ulf_PGLASS.AddVariables(fluid_model_part)   
else:
    raise "solver type not supported: "


#introducing input file name
input_file_name = problem_settings.problem_name

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


####################################################FLUID ONLY##########################################
##remove_and_save_wall_process=RemoveAndSaveWallNodesProcess()
##mark_fluid_process = MarkFluidProcess(fluid_model_part)
##
##mark_fluid_process.Execute()
##print (fluid_model_part)
##remove_and_save_wall_process.RemoveAndSave(fluid_model_part, wall_model_part)
##print (fluid_model_part)
##print (wall_model_part)

########################################################################################################

compute_reactions=problem_settings.compute_reactions
##adding dofs
if(SolverType == "Incompressible_Modified_FracStep"):
    ulf_PGLASS.AddDofs(fluid_model_part, compute_reactions)
else:
    raise "solver type not supported: "

#setting the limits of the bounding box
box_corner1 = Vector(3); 
box_corner1[0]=problem_settings.bounding_box_corner1_x; box_corner1[1]=problem_settings.bounding_box_corner1_y; box_corner1[2]=problem_settings.bounding_box_corner1_z;
box_corner2 = Vector(3); 
box_corner2[0]=problem_settings.bounding_box_corner2_x; box_corner2[1]=problem_settings.bounding_box_corner2_y; box_corner2[2]=problem_settings.bounding_box_corner2_z;

#here we write the convergence data..,
outstring2 = "convergence_info.txt"
outputfile1 = open(outstring2, 'w')

add_nodes=problem_settings.adaptive_refinement
bulk_modulus=problem_settings.bulk_modulus
density=problem_settings.density
FSI=problem_settings.FSI
#creating the solvers
#fluid solver
##check to ensure that no node has zero density or pressure

add_nodes=problem_settings.adaptive_refinement
bulk_modulus=problem_settings.bulk_modulus
density=problem_settings.density
viscosity=problem_settings.viscosity
blow_pressure=problem_settings.blow_pressure


for node in fluid_model_part.Nodes:
	node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)
	node.SetSolutionStepValue(DENSITY,0, density)	
	node.SetSolutionStepValue(VISCOSITY,0, viscosity)
	node.SetSolutionStepValue(BODY_FORCE_Y,0, -10.00000)	
	if (node.GetSolutionStepValue(FLAG_VARIABLE)==5.0):
            print ("SETTING FIXED WALL")
            node.SetSolutionStepValue(FIXED_WALL,0, 1.0)
            node.SetSolutionStepValue(FLAG_VARIABLE,0, 0.0)

for node in fluid_model_part.Nodes:
    if (node.GetSolutionStepValue(FLAG_VARIABLE)==2.0):
            node.SetSolutionStepValue(SYMMETRY_CUT,0, 1.0)
            node.SetSolutionStepValue(FLAG_VARIABLE,0, 0.0)
    #if (node.GetSolutionStepValue(FLAG_VARIABLE)==3.0):
            #node.SetSolutionStepValue(SYMMETRY_CUT,0, 2.0)
            #node.SetSolutionStepValue(FLAG_VARIABLE,0, 0.0)

is_fsi_interf=0.0
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise 'node with zero VISCOSITY found'
    if(FSI==1):
	is_fsi_interf+=node.GetSolutionStepValue(IS_INTERFACE)
	
if (SolverType == "Incompressible_Modified_FracStep" and FSI==1):
    if (is_fsi_interf==0):
	raise 'For running FSI using the Modified Frac Step Solver you must prescribe IS_INTERFACE flag at the surface/outer contour of your structure'

if(SolverType == "Incompressible_Modified_FracStep"):    
    #solver = ulf_PGLASS.ULF_FSISolver(fluid_model_part, box_corner1, box_corner2, domain_size, add_nodes, bulk_modulus, density)
    solver = ulf_PGLASS.ULF_FSISolver(fluid_model_part, structure_model_part, combined_model_part, compute_reactions, box_corner1, box_corner2, domain_size, add_nodes, blow_pressure)

    solver.alpha_shape = problem_settings.alpha_shape;
    solver.echo_level = 2;    

    solver.Initialize()
    

print "fluid solver created"


#settings to be changed
Dt = problem_settings.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01
final_time = problem_settings.max_time
output_step = problem_settings.output_step
safety_factor = 0.5 #you should put a safety factor ;-)!!!

next_output_time = output_step

time = 0.0
step = 0


outstring2 = "Acceleration.txt"
outputfile = open(outstring2, 'w')

outstring3 = "Velocity2.txt"
outputfile3 = open(outstring3, 'w')


#GlassBlowAuxProcess = GlassBlowAuxProcess(fluid_model_part);


#for node in fluid_model_part.Nodes:
	#node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)
	#node.SetSolutionStepValue(DENSITY,0, density)
	##node.SetSolutionStepValue(VISCOSITY,0, 0.2)
	#node.SetSolutionStepValue(VISCOSITY,0, 0.6)
	##node.SetSolutionStepValue(VISCOSITY,0, 4.0)
	##node.SetSolutionStepValue(VISCOSITY,0, 0.0001)
	#node.SetSolutionStepValue(BODY_FORCE_Y,0, -10.00000)	
	#if (node.GetSolutionStepValue(FLAG_VARIABLE)==5.0):
            #print ("SETTING FIXED WALL")
            #node.SetSolutionStepValue(FIXED_WALL,0, 1.0)
            #node.SetSolutionStepValue(FLAG_VARIABLE,0, 0.0)


for node in fluid_model_part.Nodes:
    
    if (node.IsFixed(DISPLACEMENT_X) and node.IsFixed(DISPLACEMENT_Y) and node.IsFixed(DISPLACEMENT_Z)):
        print "FIXED WALL"
        node.SetSolutionStepValue(FIXED_WALL,0, 1)   
    if (node.GetSolutionStepValue(FLAG_VARIABLE)==1):
        node.SetSolutionStepValue(IS_FREE_SURFACE,0, 1)
    if (node.GetSolutionStepValue(SYMMETRY_CUT)==1):
        node.Free(DISPLACEMENT_X)
        node.Free(DISPLACEMENT_Y)
        node.SetSolutionStepValue(DISPLACEMENT_Z,0, 0)
        node.Fix(DISPLACEMENT_Z)




#Hfinder  = FindNodalHProcess(fluid_model_part);
#Hfinder.Execute();   

#add_wall_process=AddWallProcess()
#counter=0

#add_wall_process.AddWall(fluid_model_part, wall_model_part)
inlet_vel = Vector(3)
dummy=LagrangianInletProcess(fluid_model_part, 0.0, inlet_vel)

while (time < final_time):
    step = step+1   
    
    
    print time
    if(step <= 3):
        new_Dt = 0.00000001;
        time = time + new_Dt*safety_factor

    #solving the fluid problem
    if(step > 3):        
        #if(time<0.002):
        #    Dt=0.001
        #elif (time<0.06):
	#    Dt=0.005
	#else:
        Dt=problem_settings.Dt 
        
        new_Dt = solver.EstimateDeltaTime(Dt, domain_size)
        #to preclude the ambiguous nodes in the part very close to the mould
        #for node in fluid_model_part.Nodes:
        #    if (node.Z>-0.035): #base nodes    and top nodes
        #        node.SetSolutionStepValue(IS_INTERFACE,0, 0.0)
                #node.Fix(DISPLACEMENT_X)
                #node.Fix(DISPLACEMENT_Z)
                #node.Fix(DISPLACEMENT_Y)
                #node.SetSolutionStepValue(FIXED_WALL,0, 1)

        #GlassBlowAuxProcess.MarkExternalSide()
        for node in fluid_model_part.Nodes:
            if (node.GetSolutionStepValue(FIXED_WALL,0)==1.0):
                node.SetSolutionStepValue(SYMMETRY_CUT,0,0.0)
       
        for node in fluid_model_part.Nodes:
            #node.SetSolutionStepValue(FLAG_VARIABLE,0, 0)
            node.SetSolutionStepValue(EXTERNAL_PRESSURE,0, 0.0)
            if (node.GetSolutionStepValue(IS_FREE_SURFACE)==1 and node.GetSolutionStepValue(IS_INTERFACE)==0):# and node.GetSolutionStepValue(SYMMETRY_CUT)==0):
                node.SetSolutionStepValue(FLAG_VARIABLE,0, 1)
                #if (time>0.2):
                node.SetSolutionStepValue(EXTERNAL_PRESSURE,0, problem_settings.blow_pressure)                
                
            if (node.GetSolutionStepValue(SYMMETRY_CUT)==1):
                "FIXING DISPL IN Z"
                node.Free(DISPLACEMENT_X)                
                node.Free(DISPLACEMENT_Y)
                node.Fix(DISPLACEMENT_Z)
           
            if (node.GetSolutionStepValue(FIXED_WALL,1)==1.0):# or node.Z>-0.015): #base nodes    and top nodes
                node.Fix(DISPLACEMENT_X)
                node.Fix(DISPLACEMENT_Z)
                node.Fix(DISPLACEMENT_Y)
                node.SetSolutionStepValue(FIXED_WALL,0, 1)

        for node in fluid_model_part.Nodes:
            if (node.GetSolutionStepValue(FIXED_WALL,0)==1.0):
                node.SetSolutionStepValue(SYMMETRY_CUT,0,0.0)

        for node in fluid_model_part.Nodes:           
            if (node.GetSolutionStepValue(SYMMETRY_CUT,0)!=0.0):                
                node.SetSolutionStepValue(EXTERNAL_PRESSURE,0, 0.0)
                node.SetSolutionStepValue(FLAG_VARIABLE,0, 0.0)  
        #density water is used to mark the inner line nodes of the symmetry cut, where external pressure must be applied
        for node in fluid_model_part.Nodes:
            if (node.GetSolutionStepValue(DENSITY_WATER)==1.0):
                #if (time>0.2):
                node.SetSolutionStepValue(EXTERNAL_PRESSURE,0, problem_settings.blow_pressure)
                node.SetSolutionStepValue(FLAG_VARIABLE,0, 1.0)  
           
                
        #new_Dt=0.001
        #impose symmetry in the middle
        for node in fluid_model_part.Nodes:
            if (node.X<0.0001 and node.X>-0.0001):
                node.Fix(DISPLACEMENT_X)
                print ("MID SURFACE - fixing displ x")
       
        time = time + new_Dt*safety_factor

        fluid_model_part.CloneTimeStep(time)
        solver.Solve(dummy)
        print "after completing the solution"

        if(time > next_output_time):
    
            file_name = input_file_name
            file_name = file_name + str(time)

            gid_io.InitializeMesh( time );
            gid_io.WriteNodeMesh((fluid_model_part).GetMesh());
            gid_io.WriteMesh((fluid_model_part).GetMesh());
            gid_io.FinalizeMesh();

            gid_io.InitializeResults(time, (fluid_model_part).GetMesh());
            #gid_io.WriteNodalResults(ACCELERATION, fluid_model_part.Nodes, time, 0);            
            #gid_io.WriteNodalResults(DISPLACEMENT, fluid_model_part.Nodes, time, 0);            
            gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, time, 0);
            #gid_io.WriteNodalResults(IS_FREE_SURFACE, fluid_model_part.Nodes, time, 0);
            #gid_io.WriteNodalResults(BODY_FORCE, fluid_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_STRUCTURE, fluid_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_INTERFACE, fluid_model_part.Nodes, time, 0);
            #gid_io.WriteNodalResults(IS_FLUID, fluid_model_part.Nodes, time, 0);
            #gid_io.WriteNodalResults(CENTER_LINE, fluid_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(PRESSURE, (fluid_model_part).Nodes, time, 0);
            #gid_io.WriteNodalResults(FIXED_WALL, (fluid_model_part).Nodes, time, 0);
            #gid_io.WriteNodalResults(NORMAL, (fluid_model_part).Nodes, time, 0);
            #gid_io.WriteNodalResults(NODAL_H, (fluid_model_part).Nodes, time, 0);
            #gid_io.WriteNodalResults(DENSITY_WATER, (fluid_model_part).Nodes, time, 0);           
            gid_io.WriteNodalResults(EXTERNAL_PRESSURE, (fluid_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(FLAG_VARIABLE, (fluid_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(SYMMETRY_CUT, (fluid_model_part).Nodes, time, 0);
            if (compute_reactions==1):
		gid_io.WriteNodalResults(REACTION, (fluid_model_part).Nodes, time, 0);
	    if (problem_settings.lagrangian_nodes_inlet==1):
                gid_io.WriteNodalResults(IS_LAGRANGIAN_INLET, (fluid_model_part).Nodes, time, 0);
                
            
            
            gid_io.Flush()
            #gid_io.CloseResultFile();
            gid_io.FinalizeResults()

            next_output_time = next_output_time  + output_step;
 
