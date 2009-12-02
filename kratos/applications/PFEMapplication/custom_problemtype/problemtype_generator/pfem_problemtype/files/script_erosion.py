#import the configuration data as read from the GiD
import pfem_antonia_var
##import edgebased_levelset_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size_pfem = pfem_antonia_var.domain_size

#read from pfem_antonia_var name and path of fluid_gid_problem
fluid_path=pfem_antonia_var.fluid_file.rpartition('/')
fluid_path=fluid_path[0]

##domain_size_fixed = edgebased_levelset_var.domain_size
##
##if (domain_size_pfem != domain_size_fixed):
##    print "Error: ------  DIFFERENT DOMAIN SIZE for the two models ------"
##else:
##    domain_size = domain_size_pfem
##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = pfem_antonia_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = pfem_antonia_var.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(fluid_path)

import edgebased_levelset_var

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_PFEMApplication = True
applications_interface.Import_MeshingApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

from KratosStructuralApplication import *
from KratosIncompressibleFluidApplication import *
from KratosPFEMApplication import *
from KratosMeshingApplication import *


## from now on the order is not anymore crucial
##################################################################
##################################################################
import math
import edgebased_levelset_solver

##import cProfile

#defining the two model part: pfem_model_part and fixed_model_part
pfem_model_part = ModelPart("PfemFluidPart");
fixed_model_part = ModelPart("FixedFluidPart");  

#importing the fluid_solver files and adding the variablesv
edgebased_levelset_solver.AddVariables(fixed_model_part)
fixed_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
##fixed_model_part.AddNodalSolutionStepVariable(DENSITY)
##fixed_model_part.AddNodalSolutionStepVariable(VISCOSITY)

#importing the structural_solver files and adding the variables
SolverType = pfem_antonia_var.SolverType
if(SolverType == "pfem_solver_ale"):
    import pfem_solver_ale_antonia
    pfem_solver_ale_antonia.AddVariables(pfem_model_part)
    pfem_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    pfem_model_part.AddNodalSolutionStepVariable(IS_VISITED)
    pfem_model_part.AddNodalSolutionStepVariable(IS_EROSIONABLE)
    pfem_model_part.AddNodalSolutionStepVariable(FRICTION_COEFFICIENT)

elif(SolverType == "monolithic_solver_lagrangian"):
    import monolithic_solver_lagrangian
    monolithic_solver_lagrangian.AddVariables(pfem_model_part)
    pfem_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    pfem_model_part.AddNodalSolutionStepVariable(IS_VISITED)
    pfem_model_part.AddNodalSolutionStepVariable(IS_EROSIONABLE)
    pfem_model_part.AddNodalSolutionStepVariable(FRICTION_COEFFICIENT)
else:
    raise "solver type not supported: options are FractionalStep - Monolithic"


#reading the models
name_pfem = pfem_antonia_var.problem_name
##print name_pfem
name_fixed = edgebased_levelset_var.problem_name
##print name_fixed

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

##READING FROM MULTIPLE .gid FILES########################

###reading the pfem model part
##gid_io = GidIO("3DDamBreak_pfem",gid_mode_flag, use_multifile, deformed_print_flag, write_conditions )
##gid_io.ReadModelPart(pfem_model_part)
##print pfem_model_part
##print "pfem model read correctly"

#reading the pfem model part
gid_io = GidIO(name_pfem,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_structure = ModelPartIO(name_pfem)
model_part_io_structure.ReadModelPart(pfem_model_part)
print "pfem model read correctly"

###reading the fixed model part with the old problem type format
##data_io = DatafileIO(name_fixed)
##data_io.ReadModelPart(fixed_model_part)
#reading the fixed model part with the new porblem type format
model_part_io_fluid = ModelPartIO(name_fixed)
model_part_io_fluid.ReadModelPart(fixed_model_part)
print "fixed model read correctly"

##NODAL CONDITIONS of the PFEM model
for node in pfem_model_part.Nodes:
    node.SetSolutionStepValue(DENSITY,0,pfem_antonia_var.Density)
    node.SetSolutionStepValue(VISCOSITY,0,pfem_antonia_var.Viscosity) 
    node.SetSolutionStepValue(BODY_FORCE_X,0,pfem_antonia_var.Gravity_X) 
    node.SetSolutionStepValue(BODY_FORCE_Y,0,pfem_antonia_var.Gravity_Y) 
    node.SetSolutionStepValue(BODY_FORCE_Z,0,pfem_antonia_var.Gravity_Z) 
    node.SetSolutionStepValue(FRICTION_COEFFICIENT,0,0.0) 


edgebased_levelset_solver.AddDofs(fixed_model_part)

##NODAL CONDITIONS of the FIXED model
##we assume here that all of the internal nodes are marked with a negative distance
##set the distance of all of the internal nodes to a small value
small_value = 0.0001
n_active = 0
for node in fixed_model_part.Nodes:
    dist = node.GetSolutionStepValue(DISTANCE)
    if(dist < 0.0):
        n_active = n_active + 1
        node.SetSolutionStepValue(DISTANCE,0,-small_value)
    else:
        node.SetSolutionStepValue(DISTANCE,0,small_value)

if(n_active == 0):
    raise "ERROR. At least one node has to be initialized with a distance lesser than 0"

###make sure that the porosity is not zero on any node (set by default to fluid only)
for node in fixed_model_part.Nodes:
    if(node.GetSolutionStepValue(POROSITY) == 0.0):
        node.SetSolutionStepValue(POROSITY,0,1.0)

##mesh_name = 0.0
##gid_io.InitializeMesh( mesh_name );
##gid_io.WriteMesh((model_part).GetMesh());
##gid_io.FinalizeMesh()

print pfem_model_part
print pfem_model_part.Properties

#setting the limits of the bounding box
box_corner1 = Vector(3);
box_corner1[0]=pfem_antonia_var.min_x;
box_corner1[1]=pfem_antonia_var.min_y;
box_corner1[2]=pfem_antonia_var.min_z;
box_corner2 = Vector(3);
box_corner2[0]=pfem_antonia_var.max_x;
box_corner2[1]=pfem_antonia_var.max_y;
box_corner2[2]=pfem_antonia_var.max_z;

#time setting
output_Dt = pfem_antonia_var.output_Dt
max_dt = pfem_antonia_var.max_dt
min_dt = pfem_antonia_var.min_dt
safety_factor = pfem_antonia_var.safety_factor
nsteps = pfem_antonia_var.nsteps

#the buffer size should be set up here after the mesh is read for the first time
pfem_model_part.SetBufferSize(2)
fixed_model_part.SetBufferSize(2)


#constructing the fluid_solver (FIXED)
body_force = Vector(3)
body_force[0] = edgebased_levelset_var.body_force_x
body_force[1] = edgebased_levelset_var.body_force_y
body_force[2] = edgebased_levelset_var.body_force_z
viscosity   = edgebased_levelset_var.viscosity
##print "*****************   VISCOSITY  *********************"
##print viscosity
density     = edgebased_levelset_var.density
##print "*****************   DENSITY  *********************"
##print density
fluid_solver = edgebased_levelset_solver.EdgeBasedLevelSetSolver(fixed_model_part,domain_size,body_force,viscosity,density)

fluid_solver.redistance_frequency = edgebased_levelset_var.redistance_frequency
fluid_solver.extrapolation_layers = edgebased_levelset_var.extrapolation_layers

fluid_solver.Initialize()
####

if(SolverType == "pfem_solver_ale"):
    #adding dofs
    pfem_solver_ale_antonia.AddDofs(pfem_model_part)
    #creating a structural_solver object
    structural_solver = pfem_solver_ale_antonia.PFEMSolver(pfem_model_part,box_corner1,box_corner2,domain_size)
    structural_solver.laplacian_form = pfem_antonia_var.laplacian_form
    structural_solver.echo_level = 0
    structural_solver.prediction_order = 1
    structural_solver.predictor_corrector = True
    structural_solver.smooth = True
    structural_solver.alpha_shape = 1.2
    structural_solver.max_press_its = 3;  
    #initializing the structural_solver
    initial_dt_pfem = 0.001*min_dt
    structural_solver.Initialize(initial_dt_pfem,output_Dt)
elif(SolverType == "monolithic_solver_lagrangian"):
    #adding dofs
    monolithic_solver_lagrangian.AddDofs(pfem_model_part)
    structural_solver = monolithic_solver_lagrangian.MonolithicSolver(pfem_model_part,domain_size,box_corner1,box_corner2)
    oss_swith = pfem_antonia_var.use_oss
    dynamic_tau = pfem_antonia_var.dynamic_tau
    pfem_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);				
    pfem_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    structural_solver.Initialize(output_Dt)


##############
ErosionModule = ErosionUtils2D();
critical_vel = Vector(3);
critical_vel[0] = 0.3;
critical_vel[1] = 0.0;
critical_vel[2] = 0.0;
critical_en = 0.1;
##############

#settings to be changed fluid FIXED
max_Dt = edgebased_levelset_var.max_time_step 
initial_Dt_ls = 0.001 * max_Dt 
final_time = edgebased_levelset_var.max_time
output_dt = edgebased_levelset_var.output_dt
safety_factor = edgebased_levelset_var.safety_factor

number_of_inital_steps = edgebased_levelset_var.number_of_inital_steps
initial_time_step = edgebased_levelset_var.initial_time_step

out = 0
time = 0.0
step = 0
next_output_time = output_dt
while(time < final_time):
    print "line49"

    if(step < number_of_inital_steps):
        Dt = initial_time_step
    else:
        Dt = fluid_solver.EstimateTimeStep(safety_factor,max_Dt)
        
    time = time + Dt
    
    fixed_model_part.CloneTimeStep(time)
    pfem_model_part.CloneTimeStep(time)

    print "******** CURRENT TIME = ",time

    if(step > 3):
        #solving the fluid problem
        fluid_solver.Solve()
##        print fixed_model_part
##        print pfem_model_part
        ErosionModule.CheckErosionableNodes(fixed_model_part, pfem_model_part, critical_vel, critical_en);

        #solving the structural problem
        structural_solver.Solve(time,gid_io)

##################################
    if(time >= next_output_time):        
        ##a change in the output name is needed!!!!
        res_name1 = str(name_fixed)
        gid_io.ChangeOutputName(res_name1)
        gid_io.InitializeMesh( time );
        gid_io.WriteMesh(fixed_model_part.GetMesh() )
        gid_io.FinalizeMesh();

        gid_io.InitializeResults( time, fixed_model_part.GetMesh() )
        gid_io.WriteNodalResults(VELOCITY,fixed_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,fixed_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISTANCE,fixed_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESS_PROJ,fixed_model_part.Nodes,time,0)
        gid_io.Flush()
        gid_io.FinalizeResults()

        
        res_name2 = str(name_pfem)
        gid_io.ChangeOutputName(res_name2)
        gid_io.InitializeMesh( time );
        gid_io.WriteNodeMesh(pfem_model_part.GetMesh());
        gid_io.WriteMesh(pfem_model_part.GetMesh() )
        gid_io.FinalizeMesh();

        gid_io.InitializeResults( time, pfem_model_part.GetMesh() )
        gid_io.WriteNodalResults(VELOCITY,pfem_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,pfem_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_FREE_SURFACE,pfem_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_BOUNDARY,pfem_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_STRUCTURE,pfem_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VISCOSITY,pfem_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DENSITY,pfem_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_FLUID,pfem_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_EROSIONABLE,pfem_model_part.Nodes,time,0)
        gid_io.Flush()  
        gid_io.FinalizeResults()

        next_output_time = time + output_dt

        out = 0

    out = out + 1
    step = step + 1
    print "step finished"

print "solution finished"        
