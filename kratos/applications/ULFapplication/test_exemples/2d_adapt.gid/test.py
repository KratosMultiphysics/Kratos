#domain size
remeshing_flag=1.0
domain_size = 2
#total simulation time
total_time = 3.0
#the max time step - it may be decreased in case it is necessary to avoid element inversion
max_delta_time = 0.003
#output time (every xxx seconds)
output_dt = 0.01
#safety factor for the delta time estimation
safety_factor = 0.5
#PATH of where the kratos library is installed
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications

project_name = '2d_adapt'

bulk_modulus = -5000.0
density = 1000.0
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
applications_interface.Import_MeshingApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)


## from now on the order is not anymore crucial
##################################################################
##################################################################
print kernel
from KratosULFApplication import *
from KratosMeshingApplication import *
from KratosStructuralApplication import *

#defining model parts
fluid_model_part = ModelPart("FluidPart");
structure_model_part = ModelPart("StructurePart");  
combined_model_part = ModelPart("CombinedPart");

fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
fluid_model_part.AddNodalSolutionStepVariable(IS_LAGRANGIAN_INLET);
fluid_model_part.AddNodalSolutionStepVariable(IS_VISITED);
fluid_model_part.AddNodalSolutionStepVariable(DISTANCE);


#adding of Variables to Model Part should be here when the "very fix container will be ready"
#importing the solver files
import ulf_fsi
ulf_fsi.AddVariables(fluid_model_part)

distance_utils=BodyDistanceCalculationUtils()

#reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(project_name,gid_mode_flag,use_multifile,deformed_print_flag,write_conditions)
gid_io.ReadModelPart(fluid_model_part)

#gid_io.ReadMesh(model_part.GetMesh())
print fluid_model_part


#the buffer size should be set up here after the mesh is read for the first time
fluid_model_part.SetBufferSize(2)

ulf_fsi.AddDofs(fluid_model_part)

#setting the limits of the bounding box
box_corner1 = Vector(3); 
box_corner1[0]= -2.5; box_corner1[1]=-2.5; box_corner1[2]=-2.5;
box_corner2 = Vector(3); 
box_corner2[0]= 3.0; box_corner2[1]=3.0; box_corner2[2]=3.0;

#creating a fluid solver object
name = project_name
solver = ulf_fsi.ULF_FSISolver(fluid_model_part, structure_model_part, combined_model_part, box_corner1, box_corner2, domain_size)
solver.alpha_shape = 1.5;
solver.echo_level = 2;


for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(BULK_MODULUS,0, bulk_modulus)
    node.SetSolutionStepValue(DENSITY,0, density)
    node.SetSolutionStepValue(VISCOSITY,0, 0.000001)
    node.SetSolutionStepValue(BODY_FORCE_Y,0, -9.8) 
        
fluid_model_part.Properties[2].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )

set_h_map_process = SetHMapProcess(fluid_model_part);

def SelectVisited(nodes):
     for node in nodes:
         if(node.GetSolutionStepValue(IS_LAGRANGIAN_INLET)==1.0):
             node.SetValue(IS_VISITED,1.0)
             node.SetSolutionStepValue(DISTANCE,0,0.0)
             #print "AAAAAAAAAAAAA"
         else:
             node.SetValue(IS_VISITED,0.0)

Dt = 0.02
nsteps = 10000
#output_Dt = 0.05
output_Dt = output_dt
min_dt = 0.00001
max_dt = max_delta_time
safety_factor = 0.5 #you should put a safety factor ;-)!!!

next_output_time = output_Dt

#initializing the solver
solver.Initialize()

SelectVisited(fluid_model_part.Nodes)


time = 0.0
step = 0

#gid_io.InitializeMesh( 0.0 );
#gid_io.WriteNodeMesh((combined_model_part).GetMesh());
#gid_io.WriteMesh((combined_model_part).GetMesh());
#gid_io.FinalizeMesh();
#gid_io.InitializeResults(0.0, (combined_model_part).GetMesh());
while (time < total_time):
    step = step+1
    
    
    
    
    print time
    if(step <= 3):
        new_Dt = 0.00000001;
        time = time + new_Dt*safety_factor

    #solving the fluid problem
    if(step > 3):
        new_Dt = solver.EstimateDeltaTime(max_dt,domain_size)
        time = time + new_Dt*safety_factor

        fluid_model_part.CloneTimeStep(time)
        structure_model_part.CloneTimeStep(time)
        combined_model_part.CloneTimeStep(time)

        print "before the solution"

        solver.Solve()

        SelectVisited(fluid_model_part.Nodes)
        distance_utils.CalculateDistances2D(fluid_model_part.Elements, DISTANCE, True)

        if (time>0.05):
            #here we want to store the min_H, max_H, min_dist, max_dist
            min_H=1000000.0
            max_H=0.0
            min_dist=0.0
            max_dist=0.0
            for node in fluid_model_part.Nodes:
                H=node.GetSolutionStepValue(NODAL_H)
                if (min_H>H):
                    min_H=H
                if (max_H<H):
                    max_H=H

            min_H=1.0*min_H
            #we set min_Dist to the min_H
            #I want to prescribe the mesh size: min_H=0.01
            #min_H=0.004
            min_Dist=0.5*min_H
            #1.0*min_H

            set_h_map_process.CalculateOptimalH(min_H, max_H)
                
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
            gid_io.WriteNodalResults(DISTANCE, (combined_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(BODY_FORCE, (combined_model_part).Nodes, time, 0);
            gid_io.FinalizeResults();

            next_output_time = next_output_time  + output_Dt;

#gid_io.FinalizeResults();

          
        


