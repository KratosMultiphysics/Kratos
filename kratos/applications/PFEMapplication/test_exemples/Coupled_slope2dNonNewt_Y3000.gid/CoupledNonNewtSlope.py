#import the configuration data as read from the GiD
import pfem_nonewtonian_coupled_var_benchmark
import time
##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = pfem_nonewtonian_coupled_var_benchmark.domain_size

#read from pfem_nonewtonian_coupled_var_benchmark name and path of fluid_gid_problem
## Acuario version*********************************
##fluid_path='/home/antonia/EXAMPLES/Erosion/WorkingErosCFixed.gid/'


##fluid_path='../Coupled_slope2dNonNewt_Y3000_fixed.gid/'
     

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking

import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_benchmarking_path)

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
import benchmarking

from KratosStructuralApplication import *
from KratosIncompressibleFluidApplication import *
from KratosPFEMApplication import *
from KratosMeshingApplication import *
from KratosExternalSolversApplication import*


## from now on the order is not anymore crucial
##################################################################
##################################################################
import math
import edgebased_levelset_solver



##import cProfile

#defining the two model part: pfem_model_part and fixed_model_part
pfem_model_part = ModelPart("PfemFluidPart");
fixed_model_part = ModelPart("FixedFluidPart");  

##Benchmarking
def NodeFinder(node_list,X,Y,Z):
   for node in node_list:
	if((node.X-X)**2 + (node.Y-Y)**2 + (node.Z-Z)**2 < .0001):
		return node

def BenchmarkCheck(time, node1):
    benchmarking.Output(time, "Time",0.1,0.0001)
##    print "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    benchmarking.Output(node1.GetSolutionStepValue(PRESSURE), "Node 1 pressure", 0.1,.0001)
    benchmarking.Output(node1.GetSolutionStepValue(VELOCITY_Y), "Node 2 velocity_y", 0.1,.0001)

#importing the fluid_solver files and adding the variablesv
edgebased_levelset_solver.AddVariables(fixed_model_part)
fixed_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
fixed_model_part.AddNodalSolutionStepVariable(POROSITY)
fixed_model_part.AddNodalSolutionStepVariable(DIAMETER)
##fixed_model_part.AddNodalSolutionStepVariable(SEEPAGE_DRAG)
##fixed_model_part.AddNodalSolutionStepVariable(VISCOSITY)

#importing the structural_solver files and adding the variables
SolverType = pfem_nonewtonian_coupled_var_benchmark.SolverType
if(SolverType == "pfem_solver_ale"):
    print "Pfem_solver_ale_not supported. Check to be using nonewtonian_asgs2d (or 3d) element"
elif(SolverType == "monolithic_solver_lagrangian"):
    import monolithic_solver_lagrangian_nonnewtonian
    monolithic_solver_lagrangian_nonnewtonian.AddVariables(pfem_model_part)
    pfem_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    pfem_model_part.AddNodalSolutionStepVariable(IS_VISITED)
    pfem_model_part.AddNodalSolutionStepVariable(IS_EROSIONABLE)
    pfem_model_part.AddNodalSolutionStepVariable(FRICTION_COEFFICIENT)
    pfem_model_part.AddNodalSolutionStepVariable(POROSITY)
    pfem_model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    pfem_model_part.AddNodalSolutionStepVariable(DIAMETER)
##    pfem_model_part.AddNodalSolutionStepVariable(SEEPAGE_DRAG)
    pfem_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE)
else:
    raise "solver type not supported: options are FractionalStep - Monolithic"


#reading the models
name_pfem = pfem_nonewtonian_coupled_var_benchmark.problem_name
name_fixed = pfem_nonewtonian_coupled_var_benchmark.fluid_file


gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

##READING FROM MULTIPLE .gid FILES########################
#reading the pfem model part
gid_io = EdgebasedGidIO(name_pfem,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_structure = ModelPartIO(name_pfem)
model_part_io_structure.ReadModelPart(pfem_model_part)
print "pfem model read correctly"
#reading the fixed model part with the new porblem type format
model_part_io_fluid = ModelPartIO(name_fixed)
model_part_io_fluid.ReadModelPart(fixed_model_part)
print "fixed model read correctly"

##NODAL CONDITIONS of the PFEM model
for node in pfem_model_part.Nodes:
    node.SetSolutionStepValue(DENSITY,0,pfem_nonewtonian_coupled_var_benchmark.Density)
    node.SetSolutionStepValue(VISCOSITY,0,pfem_nonewtonian_coupled_var_benchmark.Viscosity)
    node.SetSolutionStepValue(POROSITY,0,pfem_nonewtonian_coupled_var_benchmark.Porosity)
    node.SetSolutionStepValue(DIAMETER,0,pfem_nonewtonian_coupled_var_benchmark.Diameter)
    node.SetSolutionStepValue(BODY_FORCE_X,0,pfem_nonewtonian_coupled_var_benchmark.Gravity_X) 
    node.SetSolutionStepValue(BODY_FORCE_Y,0,pfem_nonewtonian_coupled_var_benchmark.Gravity_Y) 
    node.SetSolutionStepValue(BODY_FORCE_Z,0,pfem_nonewtonian_coupled_var_benchmark.Gravity_Z) 
    node.SetSolutionStepValue(YIELD_STRESS,0,pfem_nonewtonian_coupled_var_benchmark.YieldStress) 


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
    if(node.GetSolutionStepValue(DIAMETER) == 0.0):
        node.SetSolutionStepValue(DIAMETER,0,1.0)

##print pfem_model_part
##print pfem_model_part.Properties

#setting the limits of the bounding box
box_corner1 = Vector(3);
box_corner1[0]=pfem_nonewtonian_coupled_var_benchmark.min_x;
box_corner1[1]=pfem_nonewtonian_coupled_var_benchmark.min_y;
box_corner1[2]=pfem_nonewtonian_coupled_var_benchmark.min_z;
box_corner2 = Vector(3);
box_corner2[0]=pfem_nonewtonian_coupled_var_benchmark.max_x;
box_corner2[1]=pfem_nonewtonian_coupled_var_benchmark.max_y;
box_corner2[2]=pfem_nonewtonian_coupled_var_benchmark.max_z;

#time setting
output_dt = pfem_nonewtonian_coupled_var_benchmark.output_Dt
max_dt = pfem_nonewtonian_coupled_var_benchmark.max_dt
min_dt = pfem_nonewtonian_coupled_var_benchmark.min_dt
safety_factor = pfem_nonewtonian_coupled_var_benchmark.safety_factor
nsteps = pfem_nonewtonian_coupled_var_benchmark.nsteps

#the buffer size should be set up here after the mesh is read for the first time
pfem_model_part.SetBufferSize(2)
fixed_model_part.SetBufferSize(2)


#constructing the fluid_solver (FIXED)
body_force = Vector(3)
body_force[0] = pfem_nonewtonian_coupled_var_benchmark.Gravity_X
body_force[1] = pfem_nonewtonian_coupled_var_benchmark.Gravity_Y
body_force[2] = pfem_nonewtonian_coupled_var_benchmark.Gravity_Z

viscosity   = edgebased_levelset_var.viscosity
##print "*****************   VISCOSITY  *********************"
##print viscosity
density     = edgebased_levelset_var.density
##print "*****************   DENSITY  *********************"
##print density
fluid_solver = edgebased_levelset_solver.EdgeBasedLevelSetSolver(fixed_model_part,domain_size,body_force,viscosity,density)

fluid_solver.redistance_frequency = edgebased_levelset_var.redistance_frequency
fluid_solver.extrapolation_layers = edgebased_levelset_var.extrapolation_layers
fluid_solver.stabdt_pressure_factor = edgebased_levelset_var.stabdt_pressure_factor
fluid_solver.stabdt_convection_factor = edgebased_levelset_var.stabdt_convection_factor
fluid_solver.tau2_factor = edgebased_levelset_var.tau2_factor
fluid_solver.edge_detection_angle = edgebased_levelset_var.edge_detection_angle
fluid_solver.assume_constant_pressure = edgebased_levelset_var.assume_constant_pressure

fluid_solver.Initialize()

if(edgebased_levelset_var.wall_law_y > 1e-10):
    fluid_solver.fluid_solver.ActivateWallResistance(edgebased_levelset_var.wall_law_y);
    
####

if(SolverType == "pfem_solver_ale"):
    #adding dofs
    print "pfem_solver_ale is being used whereas with erosion module only monolithic solver lagrangian is allowed. Check also to be using use asgs2d!!"
elif(SolverType == "monolithic_solver_lagrangian"):
    #adding dofs
    monolithic_solver_lagrangian_nonnewtonian.AddDofs(pfem_model_part)
    structural_solver = monolithic_solver_lagrangian_nonnewtonian.MonolithicSolver(pfem_model_part,domain_size,box_corner1,box_corner2)
    oss_swith = pfem_nonewtonian_coupled_var_benchmark.use_oss
    dynamic_tau = pfem_nonewtonian_coupled_var_benchmark.dynamic_tau
    pfem_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);				
    pfem_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    structural_solver.linear_solver = SuperLUSolver()
##    structural_solver.conv_criteria = UPCriteria(1e-9,1e-9,1e-9,1e-9)
    structural_solver.Initialize(output_dt)



if (domain_size == 3):
    ProjectionModule = MeshTransfer3D();
    ErosionModule = ErosionUtils3D();
else:
    ProjectionModule = MeshTransfer2D();
    ErosionModule = ErosionUtils2D();
PfemUtils = PfemUtils();
DragUtils = DragUtils();

############## To be inserted in the problem type########
critical_vel = Vector(3);
critical_vel[0] = 0.06;
critical_vel[1] = 0.0;
critical_vel[2] = 0.0;
critical_en = 0.000001;
##############

ProjectionModule.DirectScalarVarInterpolation(pfem_model_part, fixed_model_part, POROSITY, POROSITY);
ProjectionModule.DirectScalarVarInterpolation(pfem_model_part, fixed_model_part, DIAMETER, DIAMETER);
for node in fixed_model_part.Nodes:
   if(node.GetSolutionStepValue(POROSITY) == 0.0):
       node.SetSolutionStepValue(POROSITY,0,1.0)
   if(node.GetSolutionStepValue(DIAMETER) == 0.0):
       node.SetSolutionStepValue(DIAMETER,0,1.0)


##bool to decide to solve or not the pfem model 
calculate_dam = True

#settings to be changed edgebased FIXED
max_Dt = edgebased_levelset_var.max_time_step
initial_Dt_ls = 0.001 * max_Dt
final_time = edgebased_levelset_var.max_time
safety_factor = edgebased_levelset_var.safety_factor
number_of_inital_steps = edgebased_levelset_var.number_of_inital_steps
initial_time_step = edgebased_levelset_var.initial_time_step
out = 0

original_max_dt = max_Dt
Time = 0.0
step = 0
next_output_time = output_dt
max_safety_factor = safety_factor

substep_number = 10
substep_counter = 0
###############################################################
back_node = NodeFinder(pfem_model_part.Nodes , 0.34001, 0.381552, 0.0)
print back_node

###############################################################
while(Time < final_time):
##    print "line49"

    if(step < number_of_inital_steps):
        max_Dt = initial_time_step
    else:
        max_Dt = original_max_dt
        #progressively increment the safety factor
        #in the steps that follow a reduction of it
        safety_factor = safety_factor * 1.2
        if(safety_factor > max_safety_factor):
            safety_factor = max_safety_factor


    Dt = fluid_solver.EstimateTimeStep(safety_factor,max_Dt)
        
    Time = Time + Dt
    fixed_model_part.CloneTimeStep(Time)

    ##let's do this later ... only each substep_number fluid solutions
    ##pfem_model_part.CloneTimeStep(Time)



    print "******** CURRENT TIME = ",Time

    if(step > 3):
        if( calculate_dam == True):
            ##(origin_model_part, destination_model_part, origin_variable_destination_variable)
            ProjectionModule.DirectScalarVarInterpolation(pfem_model_part, fixed_model_part, POROSITY, POROSITY);
            ProjectionModule.DirectScalarVarInterpolation(pfem_model_part, fixed_model_part, DIAMETER, DIAMETER);

        for node in fixed_model_part.Nodes:
            if(node.GetSolutionStepValue(POROSITY) == 0.0):
                node.SetSolutionStepValue(POROSITY,0,1.0)
            if(node.GetSolutionStepValue(DIAMETER) == 0.0):
                node.SetSolutionStepValue(DIAMETER,0,1.0)
        
        #solving the fluid problem ...........................................................................................
        print "starting solving fluid edgebased~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        fluid_solver.Solve();
        print "finished solving fluid edgebased"
        
        
        
        print "start checking the time stepping"        
        check_dt = fluid_solver.EstimateTimeStep(0.95,max_Dt)

        if(check_dt < Dt):
            print "***********************************************************"
            print "***********************************************************"
            print "***********************************************************"
            print "            *** REDUCING THE TIME STEP ***"
            print "***********************************************************"
            print "***********************************************************"
            print "***********************************************************"
                       
            #we found a velocity too large! we need to reduce the time step
            fluid_solver.fluid_solver.ReduceTimeStep(fixed_model_part,Time) ##this is to set the database to the value at the beginning of the step

            safety_factor *= edgebased_levelset_var.reduction_on_failure
            reduced_dt = fluid_solver.EstimateTimeStep(safety_factor,max_Dt)
            
            print "time before reduction= ",Time
            Time = Time - Dt + reduced_dt
            print "reduced time = ",Time
            print "Dt = ",Dt
            print "reduced_dt = ",reduced_dt

            fluid_solver.fluid_solver.ReduceTimeStep(fixed_model_part,Time) ##this is to set the database to the value at the beginning of the step
            fluid_solver.Solve()
        print "finished checking the time stepping"
        print fixed_model_part

            
        
        if(substep_counter == substep_number):
            substep_counter = 0
            
            pfem_model_part.CloneTimeStep(Time)
                       
            print pfem_model_part
            #solving the structural problem ...........................................................................................
            print "starting solving structure~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            PfemUtils.CalculateNodalMass(pfem_model_part,domain_size)
##COMMENTING OUT provisionally BEGIN
##            print "starting EROSION"
##            structural_solver.neigh_finder.Execute(); 
##            print "founded neighbours"
##            calculate_dam = ErosionModule.CheckErosionableNodes(fixed_model_part, pfem_model_part, critical_vel, critical_en, viscosity, density);
##            print "finished EROSION"
##COMMENTING OUT provisionally END
            
            if(calculate_dam == True):
                ##Projecting the water pressure
                ProjectionModule.DirectScalarVarInterpolation(fixed_model_part, pfem_model_part, PRESSURE, WATER_PRESSURE);

##COMMENTING OUT provisionally BEGIN             
                #Adding Darcy non linear term to the external forces
##                DragUtils.CalculateFluidDrag(fixed_model_part, SEEPAGE_DRAG, viscosity)
##                ProjectionModule.DirectVectorialVarInterpolation(fixed_model_part, pfem_model_part, SEEPAGE_DRAG, SEEPAGE_DRAG);
##                DragUtils.AddDrag(pfem_model_part, SEEPAGE_DRAG, BODY_FORCE, body_force)
##COMMENTING OUT provisionally END
                
                #PFEM steps
                structural_solver.Remesh();
                print "Remesh in done!"
                (structural_solver.solver).Solve();
                BenchmarkCheck(Time, back_node)
                
                print "Solve in done!"
                (structural_solver.solver).Clear();
        
                
                print "finished solving structure"
        ##        ErosionModule.SetErosionableNodes(pfem_model_part);
        ##        if (Time>0.25):
        ##        ProjectionModule.DirectScalarVarInterpolation(pfem_model_part, fixed_model_part, POROSITY);
        ##        for node in fixed_model_part.Nodes:
        ##            if( node.GetSolutionStepValue(POROSITY) == 0.5):
        ##                print node.Id

    	        if (Time > 9.0):
                    substep_number = 2
        substep_counter += 1

##################################
    if(Time >= next_output_time):        
        ##a change in the output name is needed!!!!
        
        res_name1 = str(name_fixed)
        gid_io.ChangeOutputName(res_name1)
        gid_io.InitializeMesh( Time );
        gid_io.WriteMesh(fixed_model_part.GetMesh() )
        gid_io.FinalizeMesh();
	gid_io.InitializeResults( Time, fixed_model_part.GetMesh() )
## 	    print "structure output ****************************************************"
   	gid_io.WriteNodalResults(VELOCITY,fixed_model_part.Nodes,Time,0)
## 	    gid_io.WriteNodalResults(CONV_VELOCITY,fixed_model_part.Nodes,Time,0)
   	gid_io.WriteNodalResults(PRESSURE,fixed_model_part.Nodes,Time,0)
   	gid_io.WriteNodalResults(DISTANCE,fixed_model_part.Nodes,Time,0)
   	gid_io.WriteNodalResults(PRESS_PROJ,fixed_model_part.Nodes,Time,0)
   	gid_io.WriteNodalResults(POROSITY,fixed_model_part.Nodes,Time,0)
   	gid_io.WriteNodalResults(DIAMETER,fixed_model_part.Nodes,Time,0)
##   	gid_io.WriteNodalResults(SEEPAGE_DRAG,fixed_model_part.Nodes,Time,0)


        gid_io.Flush()
        gid_io.FinalizeResults()

        
        res_name2 = str(name_pfem)
        gid_io.ChangeOutputName(res_name2)
        gid_io.InitializeMesh( Time );
        gid_io.WriteNodeMesh(pfem_model_part.GetMesh());
        gid_io.WriteMesh(pfem_model_part.GetMesh() )
        gid_io.FinalizeMesh();

        gid_io.InitializeResults( Time, pfem_model_part.GetMesh() )
##        print "fluid output ****************************************************"
        gid_io.WriteNodalResults(VELOCITY,pfem_model_part.Nodes,Time,0)
        gid_io.WriteNodalResults(PRESSURE,pfem_model_part.Nodes,Time,0)
        gid_io.WriteNodalResults(IS_FREE_SURFACE,pfem_model_part.Nodes,Time,0)
        gid_io.WriteNodalResults(IS_BOUNDARY,pfem_model_part.Nodes,Time,0)
        gid_io.WriteNodalResults(IS_STRUCTURE,pfem_model_part.Nodes,Time,0)
        gid_io.WriteNodalResults(VISCOSITY,pfem_model_part.Nodes,Time,0)
        gid_io.WriteNodalResults(DENSITY,pfem_model_part.Nodes,Time,0)
        gid_io.WriteNodalResults(BODY_FORCE,pfem_model_part.Nodes,Time,0)
##        gid_io.WriteNodalResults(FRICTION_COEFFICIENT,pfem_model_part.Nodes,Time,0)
##        gid_io.WriteNodalResults(IS_EROSIONABLE,pfem_model_part.Nodes,Time,0)
        gid_io.WriteNodalResults(POROSITY,pfem_model_part.Nodes,Time,0)
        gid_io.WriteNodalResults(DIAMETER,pfem_model_part.Nodes,Time,0)
##        gid_io.WriteNodalResults(SEEPAGE_DRAG,pfem_model_part.Nodes,Time,0)
        gid_io.WriteNodalResults(WATER_PRESSURE,pfem_model_part.Nodes,Time,0)
        gid_io.PrintOnGaussPoints(TEMPERATURE,pfem_model_part,Time)
        gid_io.PrintOnGaussPoints(AUX_INDEX,pfem_model_part,Time)

        gid_io.Flush()  
        gid_io.FinalizeResults()

        next_output_time = Time + output_dt

        out = 0

    out = out + 1
    step = step + 1
    print "step finished"

print "solution finished"        
