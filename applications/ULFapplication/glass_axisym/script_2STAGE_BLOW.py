import problem_settings
#import distance_to_wall
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
#from KratosMultiphysics.PFEMApplication import PfemUtils
from KratosMultiphysics.StructuralApplication import *

#import NistParameters
#THERMAL
#from KratosMultiphysics.ConvectionDiffusionApplication import *
#from KratosMultiphysics.ExternalSolversApplication import *

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");
structure_model_part = ModelPart("StructurePart");  
combined_model_part = ModelPart("CombinedPart");
wall_model_part = ModelPart("WallPart");

SolverType=problem_settings.SolverType

#############################################
##importing the solvers needed
if(SolverType == "Quasi_Inc_Constant_Pressure"):
  import ulf_PGLASS
  ulf_PGLASS.AddVariables(fluid_model_part)
else:
    raise "solver type not supported for GLASS simulation:can be only quasi_inc_constant_pres at the moment"

fluid_model_part.AddNodalSolutionStepVariable(TEMPERATURE)

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
compute_reactions=False

####################################################FLUID ONLY##########################################
remove_and_save_wall_process=RemoveAndSaveWallNodesProcess()
mark_fluid_process = MarkFluidProcess(fluid_model_part)

mark_fluid_process.Execute()
print (fluid_model_part)
remove_and_save_wall_process.RemoveAndSave(fluid_model_part, wall_model_part)
print (fluid_model_part)
print (wall_model_part)

########################################################################################################


#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

ulf_PGLASS.AddDofs(fluid_model_part, compute_reactions)

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
viscosity=problem_settings.viscosity
blow_pressure=problem_settings.blow_pressure

print ("LALALA")
#creating the solvers
#fluid solver
##check to ensure that no node has zero density or pressure

for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)   
    node.SetSolutionStepValue(DENSITY,0, density)
    node.SetSolutionStepValue(VISCOSITY,0, viscosity)
    node.SetSolutionStepValue(BODY_FORCE_Y,0, -9.8)   

n_flag=0
for node in fluid_model_part.Nodes:
  n_flag=n_flag+node.GetSolutionStepValue(FLAG_VARIABLE)
  
if (n_flag==0):
   print ("No nodes are marked with FLAG_VARIABLE=1, which is the inblow")
   

solver = ulf_PGLASS.ULF_FSISolver(fluid_model_part, structure_model_part, combined_model_part, compute_reactions, box_corner1, box_corner2, domain_size, add_nodes, blow_pressure)
   


solver.alpha_shape = problem_settings.alpha_shape;
solver.echo_level = 2;
solver.Initialize()

print ("fluid solver created")

is_fsi_interf=0.0

for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print(("node ",node.Id," has zero density!"))
        raise ('node with zero density found')
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print(("node ",node.Id," has zero viscosity!"))
        raise ('node with zero VISCOSITY found')
    if(node.GetSolutionStepValue(BULK_MODULUS) == 0.0):
        print(("node ",node.Id," has zero bulk modulus!"))
        raise ('node with zero BULK_MODULUS found')
    

#settings to be changed
Dt = problem_settings.Dt 
full_Dt = Dt 
initial_Dt = 0.001;# * full_Dt #0.05 #0.01
final_time = problem_settings.max_time
output_step = problem_settings.output_step
safety_factor = 0.5 #you should put a safety factor ;-)!!!

next_output_time = output_step

time = 0.0
step = 0

inlet_vel = Vector(3)

dummy=LagrangianInletProcess(fluid_model_part, 0.0, inlet_vel)

out_file = open("distance_to_wall.out", 'w')

#first blow temperature
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(TEMPERATURE,0,970.0)
    node.SetSolutionStepValue(VISCOSITY, 0, 4.0)	

#NistParameters.CalculateViscosity(fluid_model_part)



for node in fluid_model_part.Nodes:
  #if (node.GetSolutionStepValue(IS_STRUCTURE)==1.0 and node.GetSolutionStepValue(IS_FLUID)!=1.0):
  #  node.SetSolutionStepValue(NODAL_H, 0, 0.0005)
  #node.SetSolutionStepValue(IS_STRUCTURE, 0, 0.0)
  node.SetSolutionStepValue(BODY_FORCE_Y, 0, -9.8)
  if (node.Y<0.005 and node.X>0.006):
    node.SetSolutionStepValue(VISCOSITY, 0, 40.0)

for node in wall_model_part.Nodes:
  node.SetSolutionStepValue(IS_STRUCTURE, 0, 1)

add_wall_process=AddWallProcess()
counter=0

while (time < final_time):
    step = step+1   
    
    
    print (time)
    if(step <= 3):
        new_Dt = 0.00000001;
        time = time + new_Dt*safety_factor

    #solving the fluid problem
    if(step > 3):
        new_Dt = solver.EstimateDeltaTime(Dt, domain_size)
        print((("DT IS "), Dt))
        time = time + new_Dt*safety_factor

        combined_model_part.CloneTimeStep(time)
        fluid_model_part.ProcessInfo=combined_model_part.ProcessInfo
        structure_model_part.ProcessInfo=combined_model_part.ProcessInfo	  
        
        #if (time>0.4 and counter==0):
        #for measured temps/visc
        if (time>1.6 and counter==0):
          #first blow temperature
##          for node in fluid_model_part.Nodes:
##            #node.SetSolutionStepValue(TEMPERATURE,0,750.0)
##            if (node.Y<0.011 and node.X>0.006):
##              node.SetSolutionStepValue(VISCOSITY, 0, 100.0)
##            #VISCOSIDAD PROPORCIONADA POR RC
##            else:
##              node.SetSolutionStepValue(VISCOSITY, 0, 40.0)
              
          for node in fluid_model_part.Nodes:
            if (node.Y>0.013): #THATS WHERE SECOND MOULD STARTS
              if (node.GetSolutionStepValue(IS_STRUCTURE)==1):
                node.SetSolutionStepValue(IS_INTERFACE,0, 1)
                node.SetSolutionStepValue(IS_STRUCTURE,0, 0)
              node.Free(DISPLACEMENT_X)
              node.Free(DISPLACEMENT_Y)
              
              
	  add_wall_process.AddWall(fluid_model_part, wall_model_part)
	  
	  for node in fluid_model_part.Nodes: 
              node.SetSolutionStepValue(VISCOSITY, 0, 40.0)
              node.SetSolutionStepValue(BODY_FORCE_Y, 0, 9.8)
	  
	  counter=1
       
        solver.Solve(dummy)
        #in the last step we write the thickness ditribution
        if (time>final_time-full_Dt):
            out_file.write( str(time)+ "\n")
            distance_to_wall.compute_distance_to_wall(combined_model_part, out_file)
                      
        print ("after completing the solution")

        if(time > next_output_time):
    
            file_name = input_file_name
            file_name = file_name + str(time)

            gid_io.InitializeMesh( time );
            gid_io.WriteNodeMesh((combined_model_part).GetMesh());
            gid_io.WriteMesh((combined_model_part).GetMesh());
            gid_io.FinalizeMesh();

            gid_io.InitializeResults(time, (combined_model_part).GetMesh());
    
            gid_io.WriteNodalResults(TEMPERATURE, combined_model_part.Nodes, time, 0);            
            gid_io.WriteNodalResults(VELOCITY, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(PRESSURE, (combined_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(EXTERNAL_PRESSURE, (combined_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(IS_INTERFACE, (combined_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(IS_STRUCTURE, (combined_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(VISCOSITY, (combined_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(IS_FLUID, (combined_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(BODY_FORCE, (combined_model_part).Nodes, time, 0);            
            if (compute_reactions==1):
                gid_io.WriteNodalResults(REACTION, (combined_model_part).Nodes, time, 0);           
            gid_io.Flush()
            #gid_io.CloseResultFile();
            gid_io.FinalizeResults()

            next_output_time = next_output_time  + output_step;
 
