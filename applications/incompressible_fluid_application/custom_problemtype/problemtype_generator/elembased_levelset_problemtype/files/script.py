import elembased_levelset_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = elembased_levelset_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = elembased_levelset_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = elembased_levelset_var.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_ConvectionDiffusionApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosIncompressibleFluidApplication import *
from KratosConvectionDiffusionApplication import *

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  

#############################################


##importing the solvers needed
import level_set_elembased_fluid_solver
level_set_elembased_fluid_solver.AddVariables(fluid_model_part)

#introducing input file name
input_file_name = elembased_levelset_var.problem_name

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

##selecting output format
if(elembased_levelset_var.print_layers == True):
    gid_io = EdgebasedGidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
else:
    gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
    
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(2)

##adding dofs
level_set_elembased_fluid_solver.AddDofs(fluid_model_part)

##we assume here that all of the internal nodes are marked with a negative distance
##set the distance of all of the internal nodes to a small value
small_value = 0.0001
n_active = 0
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,elembased_levelset_var.viscosity)
    node.SetSolutionStepValue(DENSITY,0,elembased_levelset_var.density)
    node.Free(PRESSURE)
    node.SetSolutionStepValue(PRESSURE,0,0.0)
    dist = node.GetSolutionStepValue(DISTANCE)
    if(dist < 0.0):
        n_active = n_active + 1
        node.SetSolutionStepValue(DISTANCE,0,-small_value)
    else:
        node.SetSolutionStepValue(DISTANCE,0,small_value)

if(n_active == 0):
    raise "ERROR. At least one node has to be initialized with a distance lesser than 0"

#make sure that the porosity is not zero on any node (set by default to fluid only)
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(POROSITY) == 0.0):
        node.SetSolutionStepValue(POROSITY,0,1.0)
    if(node.GetSolutionStepValue(DIAMETER) == 0.0):
        node.SetSolutionStepValue(DIAMETER,0,1.0)
        
#constructing the solver
body_force = Vector(3)
body_force[0] = elembased_levelset_var.body_force_x
body_force[1] = elembased_levelset_var.body_force_y
body_force[2] = elembased_levelset_var.body_force_z
viscosity   = elembased_levelset_var.viscosity
density     = elembased_levelset_var.density
fluid_solver = level_set_elembased_fluid_solver.ElemBasedLevelSetSolver(fluid_model_part,domain_size,body_force)

fluid_solver.redistance_frequency = elembased_levelset_var.redistance_frequency
fluid_solver.extrapolation_layers = elembased_levelset_var.extrapolation_layers

fluid_solver.Initialize()
####


print "fluid solver created"

#settings to be changed
Dt = elembased_levelset_var.time_step  
final_time = elembased_levelset_var.max_time
output_dt = elembased_levelset_var.output_dt
coef = elembased_levelset_var.delta_time_coefficient

#number_of_inital_steps = elembased_levelset_var.number_of_inital_steps
#initial_time_step = elembased_levelset_var.initial_time_step
out = 0


###mesh to be printed
if(elembased_levelset_var.print_layers == False):
    mesh_name = 0.0
    gid_io.InitializeMesh( mesh_name)
    gid_io.WriteMesh( fluid_model_part.GetMesh() )
    gid_io.FinalizeMesh()
    gid_io.Flush()

    gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh());

time_old_print = 0.0    
time = 0.0
step = 0

next_output_time = output_dt
while(time < final_time):

        
    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    print "******** CURRENT TIME = ",time

    if(step >= 3):
        ##Calculate Dt when a jump in velocity is reached
        Dt_old = elembased_levelset_var.time_step
        Dt_new = fluid_solver.CalculateDelta_t(Dt)
        if(Dt_old >= coef * Dt_new):
            Dt = coef * Dt_new

        fluid_solver.Solve()

    time_to_print = time - time_old_print
#    if(time >= next_output_time):
    if(time_to_print >= elembased_levelset_var.output_dt ):
        if(elembased_levelset_var.print_layers == True):
            #writing mesh 
            gid_io.InitializeMesh( time );
            gid_io.WriteMesh((fluid_model_part).GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(time, (fluid_model_part).GetMesh());
            
        gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISTANCE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VISCOSITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DENSITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(NORMAL,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(POROSITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DIAMETER,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_STRUCTURE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(CONVECTION_VELOCITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(BODY_FORCE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(AUX_INDEX,fluid_model_part.Nodes,time,0)
        gid_io.Flush()

        if(elembased_levelset_var.print_layers == True):
            gid_io.FinalizeResults()
	time_old_print = time
#        next_output_time = time + output_dt


 #       out = 0

 #   out = out + 1
    step = step + 1
      
if(elembased_levelset_var.print_layers == False):
    gid_io.FinalizeResults()    
        

