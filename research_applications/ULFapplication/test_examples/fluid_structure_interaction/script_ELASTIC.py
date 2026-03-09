import fluid_ulf_var
##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = fluid_ulf_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
#importing Kratos main library
#import sys
import math
import sys
sys.path.append(fluid_ulf_var.kratos_path)
from KratosMultiphysics import *
#from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
#from KratosMultiphysics.PFEMApplication import PfemUtils
#STR APPL U NEED FOR COMPUTING DISTANCES BETWEEN VOLUME AND SURFACE MESHES
from KratosMultiphysics.StructuralMechanicsApplication import *



#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");
structure_model_part = ModelPart("StructurePart");  
combined_model_part = ModelPart("CombinedPart");

SolverType=fluid_ulf_var.SolverType
if (SolverType=="Incompressible_Modified_FracStep" or SolverType=="FracStep"):
    fluid_only_model_part = ModelPart("FluidOnlyPart");

#############################################
##importing the solvers needed
if(SolverType == "Incompressible_Modified_FracStep"):
    import ulf_frac
    ulf_frac.AddVariables(fluid_model_part)   
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
    raise ("solver type not supported:")

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

compute_reactions=fluid_ulf_var.compute_reactions
##adding dofs
if(SolverType == "Incompressible_Modified_FracStep"):
    ulf_frac.AddDofs(fluid_model_part, compute_reactions)
elif(SolverType == "FracStep"):
    ulf_frac.AddDofs(fluid_model_part, compute_reactions)
elif(SolverType == "Quasi_Inc_Constant_Pressure"):
    ulf_fsi.AddDofs(fluid_model_part, compute_reactions)
elif(SolverType == "Quasi_Inc_Linear_Pressure"):
    ulf_fsi_inc.AddDofs(fluid_model_part, compute_reactions)






#setting the limits of the bounding box
box_corner1 = Vector(3); 
box_corner1[0]=fluid_ulf_var.bounding_box_corner1_x; box_corner1[1]=fluid_ulf_var.bounding_box_corner1_y; box_corner1[2]=fluid_ulf_var.bounding_box_corner1_z;
box_corner2 = Vector(3); 
box_corner2[0]=fluid_ulf_var.bounding_box_corner2_x; box_corner2[1]=fluid_ulf_var.bounding_box_corner2_y; box_corner2[2]=fluid_ulf_var.bounding_box_corner2_z;

#here we write the convergence data..,
outstring2 = "convergence_info.txt"
outputfile1 = open(outstring2, 'w')

add_nodes=fluid_ulf_var.adaptive_refinement
bulk_modulus=fluid_ulf_var.bulk_modulus
density=fluid_ulf_var.density
viscosity=fluid_ulf_var.viscosity
body_force=Vector(3)

if (fluid_ulf_var.with_gravity==1):
    if (domain_size==2):
      body_force[1]=-9.8
    elif (domain_szie==3):
      body_force[2]=-9.8

FSI=fluid_ulf_var.FSI

#creating the solvers
#fluid solver
##check to ensure that no node has zero density or pressure
is_fsi_interf=0.0

for node in fluid_model_part.Nodes:
    if (node.GetSolutionStepValue(IS_STRUCTURE)==1.0 and node.IsFixed(DISPLACEMENT_X)==False and node.IsFixed(DISPLACEMENT_Y)==False and node.IsFixed(DISPLACEMENT_Z)==False):
      node.SetSolutionStepValue(IS_INTERFACE,0,1.0)


for node in fluid_model_part.Nodes:
    if(FSI==1):
      is_fsi_interf+=node.GetSolutionStepValue(IS_INTERFACE)      
if (SolverType == "Incompressible_Modified_FracStep" and FSI==1):
    if (is_fsi_interf==0):
      raise ('For running FSI using the Modified Frac Step Solver you must prescribe IS_INTERFACE flag at the surface/outer contour of your structure')

if(SolverType == "Incompressible_Modified_FracStep"):    
    solver = ulf_frac.ULF_FSISolver(outputfile1, fluid_only_model_part, fluid_model_part, structure_model_part, combined_model_part, FSI, compute_reactions, box_corner1, box_corner2, domain_size, add_nodes, bulk_modulus, density)
    solver.alpha_shape = fluid_ulf_var.alpha_shape;
    solver.echo_level = 2;
    
    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)
        node.SetSolutionStepValue(DENSITY,0, density)   
        node.SetSolutionStepValue(VISCOSITY,0, viscosity)   
        node.SetSolutionStepValue(BODY_FORCE,0, body_force)
    solver.Initialize()


print ("fluid solver created")

for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print ("node ",node.Id," has zero density!")
        raise ('node with zero density found')
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print ("node ",node.Id," has zero viscosity!")
        raise ('node with zero VISCOSITY found')    

if (fluid_ulf_var.FSI==1):
    if (fluid_ulf_var.domain_size==2):
      prop = fluid_model_part.Properties[1]
      mat =LinearElasticPlaneStrain2DLaw()#LinearElasticPlaneStrain2DLaw()#HyperElasticPlaneStrain2DLaw();
      prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());	
    elif (fluid_ulf_var.domain_size==3):
      fluid_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
    else:
      raise "Domain size error. It should be 2D or 3D"

print ()

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

dummy=LagrangianInletProcess(fluid_model_part, 0.0, inlet_vel)


while (time < final_time):
    step = step+1   
    
    
    print (time)
    if(step <= 3):
        new_Dt = 0.00000001;
        time = time + new_Dt*safety_factor

    #solving the fluid problem
    if(step > 3):
        new_Dt = solver.EstimateDeltaTime(Dt, domain_size)
        time = time + new_Dt*safety_factor
        fluid_model_part.CloneTimeStep(time)  
        combined_model_part.ProcessInfo=fluid_model_part.ProcessInfo
        structure_model_part.ProcessInfo=fluid_model_part.ProcessInfo
        
        solver.Solve(dummy)
               
        print ("after completing the solution")

        if(time > next_output_time):
    
            file_name = input_file_name
            file_name = file_name + str(time)

            gid_io.InitializeMesh( time );
            gid_io.WriteNodeMesh((combined_model_part).GetMesh());
            gid_io.WriteMesh((combined_model_part).GetMesh());
            gid_io.FinalizeMesh();

            gid_io.InitializeResults(time, (combined_model_part).GetMesh());
            gid_io.WriteNodalResults(IS_STRUCTURE, combined_model_part.Nodes, time, 0);            
            gid_io.WriteNodalResults(DISPLACEMENT, combined_model_part.Nodes, time, 0);            
            gid_io.WriteNodalResults(VELOCITY, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(PRESSURE, (combined_model_part).Nodes, time, 0);
                
            
            
            gid_io.Flush()
            #gid_io.CloseResultFile();
            gid_io.FinalizeResults()

            next_output_time = next_output_time  + output_step;
 
