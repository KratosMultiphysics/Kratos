#import the configuration data as read from the GiD
import Kratos_Structural_Application_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = Kratos_Structural_Application_var.domain_size

##################################################################
##################################################################
##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path


kratos_path = '../../../../' ##kratos_root/
kratos_benchmarking_path = '../../../../benchmarking/' ##kratos_root/benchmarking
import sys
sys.path.append(kratos_path)

#importing Kratos main library
from KratosMultiphysics import *

from KratosMultiphysics.StructuralApplication import *
#loading meshing application
sys.path.append(kratos_applications_path + 'meshing_application/python_scripts') 
from KratosMultiphysics.MeshingApplication import *
meshing_application = KratosMeshingApplication()
kernel.AddApplication(meshing_application)
applications_interface.ImportApplications(kernel, kratos_applications_path)


## from now on the order is not anymore crucial
#############################/data/ArchivosKratos/pendulo/Pendulo_Kratos.gid#####################################
##################################################################

#defining a model part
model_part = ModelPart("FluidPart");  
model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
#model_part.AddNodalSolutionStepVariable(FORCE);
#model_part.AddNodalSolutionStepVariable(NODAL_DAMAGE);
#model_part.AddNodalSolutionStepVariable(NODAL_AREA);
#model_part.AddNodalSolutionStepVariable(NODAL_VALUES);
#model_part.AddNodalSolutionStepVariable(SPLIT_NODAL);
#model_part.AddNodalSolutionStepVariable(IS_DUPLICATED);
#model_part.AddNodalSolutionStepVariable(NODAL_STRAIN);
#model_part.AddNodalSolutionStepVariable(NODAL_STRESS);
#model_part.AddNodalSolutionStepVariable(DAMAGE);
#model_part.AddNodalSolutionStepVariable(DENSITY);
#model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
#model_part.AddNodalSolutionStepVariable(REACTION);
#model_part.AddNodalSolutionStepVariable(BODY_FORCE);



#reading a model
name = Kratos_Structural_Application_var.problem_name

gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(name,gid_mode_flag, use_multifile, deformed_print_flag, write_conditions)
model_part_io = ModelPartIO(name)
model_part_io.ReadModelPart(model_part)



##find neighbours if required
number_of_avg_elems = 10
number_of_avg_nodes = 10
nodal_neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
nodal_neighbour_search.Execute()
neighbour_calculator = FindElementalNeighboursProcess(model_part,2,10);
neighbour_calculator.Execute()



#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3) # me permite acceder a timepos hasta 3 pasos atras



model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
#model_part.Properties[2].SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
print "Linear elastic model selected"

for elem in model_part.Elements:
	  elem.SetValue(REFINEMENT_LEVEL,0)


model_part.Elements[11].SetValue(REFINEMENT_LEVEL,1)
model_part.Elements[12].SetValue(REFINEMENT_LEVEL,1)
model_part.Elements[13].SetValue(REFINEMENT_LEVEL,1)
model_part.Elements[47].SetValue(REFINEMENT_LEVEL,1)
model_part.Elements[48].SetValue(REFINEMENT_LEVEL,1)


Refine = LocalRefineTriangleMesh(model_part) 


Dt      = 0.0001
maxtime = 3
nsteps  = 10




for step in range(0,nsteps):


    print "TIME STEP = ", step
    model_part.ProcessInfo[TIME_STEPS] = step
    
    if(step >= 3):
      
        #refinando los elementos seleccionados
        refine_on_reference = True
        interpolate_internal_variables = False
        Refine.LocalRefineMesh(refine_on_reference,interpolate_internal_variables) 

        #recalculating neighbours 
        nodal_neighbour_search.Execute()


        gid_io.InitializeResults(step,(model_part).GetMesh())
        gid_io.InitializeMesh( step );
        gid_io.WriteNodeMesh((model_part).GetMesh());
        gid_io.WriteMesh((model_part).GetMesh());
        gid_io.FinalizeMesh();

	 
        gid_io.Flush() 
        gid_io.FinalizeResults()

print "COMPLETED ANALYSIS" 
    
        

