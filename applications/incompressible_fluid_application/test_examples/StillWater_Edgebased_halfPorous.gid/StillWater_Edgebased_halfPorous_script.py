import edgebased_levelset_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = edgebased_levelset_var.domain_size

import math
##import cProfile
##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_path = '../../../..'
kratos_benchmarking_path = '../../../../benchmarking'
import sys
sys.path.append(kratos_path)
sys.path.append(kratos_benchmarking_path)
import benchmarking


## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *

##################################################################
##################################################################
def NodeFinder(node_list,X,Y,Z):
   for node in node_list:
	if((node.X-X)**2 + (node.Y-Y)**2 + (node.Z-Z)**2 < .000001):
		return node


def BenchmarkCheck(time, node1):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(PRESSURE), "Node 1 pressure", 0.05)
    benchmarking.Output(node1.GetSolutionStepValue(DISTANCE), "Node 1 distance", 0.05)

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  

##importing the solvers needed
import edgebased_levelset_solver
edgebased_levelset_solver.AddVariables(fluid_model_part)

#introducing input file name
input_file_name = edgebased_levelset_var.problem_name

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions

##selecting output format
if(edgebased_levelset_var.print_layers == True):
    gid_io = EdgebasedGidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
else:
    gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
    
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(2)

##adding dofs
edgebased_levelset_solver.AddDofs(fluid_model_part)

##we assume here that all of the internal nodes are marked with a negative distance
##set the distance of all of the internal nodes to a small value
small_value = 0.0001
n_active = 0
for node in fluid_model_part.Nodes:
    dist = node.GetSolutionStepValue(DISTANCE)
    node.SetSolutionStepValue(DISTANCE,0,node.Y - 5.0)
    if(node.X >= 5.0):
        node.SetSolutionStepValue(POROSITY,0, 0.5)
        node.SetSolutionStepValue(DIAMETER,0, 0.35)
       

#make sure that the porosity is not zero on any node (set by default to fluid only)
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(POROSITY) == 0.0):
        node.SetSolutionStepValue(POROSITY,0,1.0)
    if(node.GetSolutionStepValue(DIAMETER) == 0.0):
        node.SetSolutionStepValue(DIAMETER,0,1.0)
        
#constructing the solver
body_force = Vector(3)
body_force[0] = edgebased_levelset_var.body_force_x
body_force[1] = edgebased_levelset_var.body_force_y
body_force[2] = edgebased_levelset_var.body_force_z
if(body_force[0] == 0.0 and body_force[1] == 0.0 and body_force[2] == 0.0):
    raise "ERROR. Body Force cannot be a ZERO VECTOR"

viscosity   = edgebased_levelset_var.viscosity
density     = edgebased_levelset_var.density
fluid_solver = edgebased_levelset_solver.EdgeBasedLevelSetSolver(fluid_model_part,domain_size,body_force,viscosity,density)
fluid_solver.redistance_frequency = edgebased_levelset_var.redistance_frequency
fluid_solver.extrapolation_layers = int(edgebased_levelset_var.extrapolation_layers)
fluid_solver.stabdt_pressure_factor = edgebased_levelset_var.stabdt_pressure_factor
fluid_solver.stabdt_convection_factor = edgebased_levelset_var.stabdt_convection_factor
fluid_solver.use_mass_correction = edgebased_levelset_var.use_mass_correction
fluid_solver.tau2_factor = edgebased_levelset_var.tau2_factor
fluid_solver.edge_detection_angle = edgebased_levelset_var.edge_detection_angle
fluid_solver.assume_constant_pressure = edgebased_levelset_var.assume_constant_pressure
fluid_solver.compute_porous_resistance_law = int(edgebased_levelset_var.compute_porous_resistance_law) ##0 = None; 1 = Ergun; 2 = Custom;
##print "compute_porous_resistance_law   ", fluid_solver.compute_porous_resistance_law
#using MKLPardisosolver ----> it has to be compiled in kratos!!
##fluid_solver.pressure_linear_solver = MKLPardisoSolver()

fluid_solver.dynamic_tau = 1.0
             
fluid_solver.Initialize()

if(edgebased_levelset_var.wall_law_y > 1e-10):
    fluid_solver.fluid_solver.ActivateWallResistance(edgebased_levelset_var.wall_law_y);
    
####


print "fluid solver created"

#settings to be changed
max_Dt = edgebased_levelset_var.max_time_step 
initial_Dt = 0.001 * max_Dt 
final_time = edgebased_levelset_var.max_time
output_dt = edgebased_levelset_var.output_dt
safety_factor = edgebased_levelset_var.safety_factor

number_of_inital_steps = edgebased_levelset_var.number_of_inital_steps
initial_time_step = edgebased_levelset_var.initial_time_step
out = 0

original_max_dt = max_Dt

###############################################################
back_node = NodeFinder(fluid_model_part.Nodes , 4.804, 2.0, 0.0)
print back_node
###############################################################

###mesh to be printed
if(edgebased_levelset_var.print_layers == False):
    mesh_name = 0.0
    gid_io.InitializeMesh( mesh_name)
    gid_io.WriteMesh( fluid_model_part.GetMesh() )
    gid_io.FinalizeMesh()
    gid_io.Flush()

    gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh());

max_safety_factor = safety_factor
    
time = 0.0
step = 0
next_output_time = output_dt
while(time < final_time):

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

        
    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    print "******** CURRENT TIME = ",time

    if(step >= 3):
        fluid_solver.Solve()
        BenchmarkCheck(time, back_node)

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
            fluid_solver.fluid_solver.ReduceTimeStep(fluid_model_part,time) ##this is to set the database to the value at the beginning of the step

            safety_factor *= edgebased_levelset_var.reduction_on_failure
            reduced_dt = fluid_solver.EstimateTimeStep(safety_factor,max_Dt)
            
            print "time before reduction= ",time
            time = time - Dt + reduced_dt
            print "reduced time = ",time
            print "Dt = ",Dt
            print "reduced_dt = ",reduced_dt

            fluid_solver.fluid_solver.ReduceTimeStep(fluid_model_part,time) ##this is to set the database to the value at the beginning of the step
           
            fluid_solver.Solve()

    if(time >= next_output_time):
        if(edgebased_levelset_var.print_layers == True):
            #writing mesh 
            gid_io.InitializeMesh( time );
            gid_io.WriteMesh((fluid_model_part).GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(time, (fluid_model_part).GetMesh());
            
        gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(POROSITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISTANCE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESS_PROJ,fluid_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(LIN_DARCY_COEF,fluid_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(NONLIN_DARCY_COEF,fluid_model_part.Nodes,time,0)
        gid_io.Flush()

        if(edgebased_levelset_var.print_layers == True):
            gid_io.FinalizeResults()

        next_output_time = time + output_dt

        out = 0

    out = out + 1
    step = step + 1
      
if(edgebased_levelset_var.print_layers == False):
    gid_io.FinalizeResults()    
        
 
