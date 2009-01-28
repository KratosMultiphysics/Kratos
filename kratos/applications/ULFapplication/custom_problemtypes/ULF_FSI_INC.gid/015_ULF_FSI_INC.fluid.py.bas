#domain size
domain_size = *GenData(DOMAIN_SIZE)
# remesh 1.0 or not remesh 0.0
remeshing_flag=*GenData(REMESHING_FLAG)
#total simulation time
total_time = *GenData(TOTAL_TIME)
#the max time step - it may be decreased in case it is necessary to avoid element inversion
max_delta_time = *GenData(MAX_DELTA_TIME)
#output time (every xxx seconds)
output_dt = *GenData(OUTPUT_DT)
#safety factor for the delta time estimation
safety_factor = *GenData(SAFETY_FACTOR)
#PATH of where the kratos library is installed
kratos_libs_path = *GenData(KRATOS_LIBS_PATH)
#PATH of where your application is installed
kratos_applications_path = *GenData(KRATOS_APPLICATIONS_PATH)
project_name = '*tcl(file tail [GiD_Info project Modelname])'

bulk_modulus = *GenData(BULK_MODULUS)
density = *GenData(DENSITY)
##################################################################
##################################################################
## ATTENTION: here the order is important

import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import **
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
from KratosR1ULFApplication import **
from KratosR1StructuralApplication import **

#defining model parts
fluid_model_part = ModelPart("FluidPart");
structure_model_part = ModelPart("StructurePart");  
combined_model_part = ModelPart("CombinedPart");

#adding of Variables to Model Part should be here when the "very fix container will be ready"
#importing the solver files
import ulf_fsi_inc
ulf_fsi_inc.AddVariables(fluid_model_part)

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

ulf_fsi_inc.AddDofs(fluid_model_part)

#setting the limits of the bounding box
box_corner1 = Vector(3); 
box_corner1[0]= *GenData(BOUNDING_BOX_LOWER_CORNER_X); box_corner1[1]=*GenData(BOUNDING_BOX_LOWER_CORNER_Y); box_corner1[2]=*GenData(BOUNDING_BOX_LOWER_CORNER_Z);
box_corner2 = Vector(3); 
box_corner2[0]= *GenData(BOUNDING_BOX_UPPER_CORNER_X); box_corner2[1]=*GenData(BOUNDING_BOX_UPPER_CORNER_Y); box_corner2[2]=*GenData(BOUNDING_BOX_UPPER_CORNER_Z);

#creating a fluid solver object
name = project_name
solver = ulf_fsi_inc.ULF_FSISolver(fluid_model_part, structure_model_part, combined_model_part, box_corner1, box_corner2, domain_size, remeshing_flag, bulk_modulus, density)
solver.alpha_shape = *GenData(ALPHA_SHAPE);
solver.echo_level = 2;

for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)
    node.SetSolutionStepValue(DENSITY,0, density) 

*loop materials
*if(strcmp(MatProp(1),"Solid")==0)
*format "%i%f"

*if(GenData(DOMAIN_SIZE,int)==2)
fluid_model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
*elseif(GenData(DOMAIN_SIZE,int)==3)
fluid_model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
*endif

*endif
*end materials


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

time = 0.0
step = 0


while (time < total_time):
    step = step+1
    
    
    
    
    print time
    if(step <= 3):
        new_Dt = 0.00000001;
	time = time + new_Dt**safety_factor

    #solving the fluid problem
    if(step > 3):
	new_Dt = solver.EstimateDeltaTime(max_dt,domain_size)
	time = time + new_Dt**safety_factor

        fluid_model_part.CloneTimeStep(time)
        structure_model_part.CloneTimeStep(time)
        combined_model_part.CloneTimeStep(time)

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
            gid_io.WriteNodalResults(FORCE, (combined_model_part).Nodes, time, 0);

            gid_io.Flush()
            #gid_io.CloseResultFile();
            gid_io.FinalizeResults()

            next_output_time = next_output_time  + output_Dt;
 

          
        


