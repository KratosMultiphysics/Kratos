#import the configuration data as read from the GiD
import pfem_nonewtonian_coupled_var
import time
##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = pfem_nonewtonian_coupled_var.domain_size

#read from pfem_nonewtonian_coupled_var name and path of fluid_gid_problem
## Acuario version*********************************
##fluid_path='/home/antonia/EXAMPLES/Erosion/WorkingErosCFixed.gid/'


fluid_path=pfem_nonewtonian_coupled_var.fluid_file.rpartition('/')
fluid_path=fluid_path[0]

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = pfem_nonewtonian_coupled_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = pfem_nonewtonian_coupled_var.kratos_path + '/applications' ##kratos_root/applications
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
applications_interface.Import_ExternalSolversAppliction = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

##from KratosStructuralApplication import *
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

#importing the fluid_solver files and adding the variablesv
edgebased_levelset_solver.AddVariables(fixed_model_part)
fixed_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
fixed_model_part.AddNodalSolutionStepVariable(POROSITY)
fixed_model_part.AddNodalSolutionStepVariable(DIAMETER)
fixed_model_part.AddNodalSolutionStepVariable(SEEPAGE_DRAG)
##fixed_model_part.AddNodalSolutionStepVariable(VISCOSITY)

#importing the structural_solver files and adding the variables
SolverType = "monolithic_solver_lagrangian"
import monolithic_solver_lagrangian_nonnewtonian
monolithic_solver_lagrangian_nonnewtonian.AddVariables(pfem_model_part)
pfem_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
pfem_model_part.AddNodalSolutionStepVariable(IS_VISITED)
pfem_model_part.AddNodalSolutionStepVariable(IS_EROSIONABLE)
pfem_model_part.AddNodalSolutionStepVariable(FRICTION_COEFFICIENT)
pfem_model_part.AddNodalSolutionStepVariable(POROSITY)
pfem_model_part.AddNodalSolutionStepVariable(NODAL_MASS)
pfem_model_part.AddNodalSolutionStepVariable(DIAMETER)
pfem_model_part.AddNodalSolutionStepVariable(SEEPAGE_DRAG)
pfem_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE)
##pfem_model_part.AddNodalSolutionStepVariable(AUX_INDEX)



#reading the models
name_pfem = pfem_nonewtonian_coupled_var.problem_name
name_fixed = pfem_nonewtonian_coupled_var.fluid_file


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
body_force = Vector(3)
body_force[0] = pfem_nonewtonian_coupled_var.Gravity_X
body_force[1] = pfem_nonewtonian_coupled_var.Gravity_Y
body_force[2] = pfem_nonewtonian_coupled_var.Gravity_Z
for node in pfem_model_part.Nodes:
    node.SetSolutionStepValue(DENSITY,0,pfem_nonewtonian_coupled_var.Density)
    node.SetSolutionStepValue(VISCOSITY,0,pfem_nonewtonian_coupled_var.Viscosity)
    node.SetSolutionStepValue(POROSITY,0,pfem_nonewtonian_coupled_var.Porosity)
    node.SetSolutionStepValue(DIAMETER,0,pfem_nonewtonian_coupled_var.Diameter)
    node.SetSolutionStepValue(BODY_FORCE_X,0,body_force[0]) 
    node.SetSolutionStepValue(BODY_FORCE_Y,0,body_force[1]) 
    node.SetSolutionStepValue(BODY_FORCE_Z,0,body_force[2]) 
    node.SetSolutionStepValue(YIELD_STRESS,0,pfem_nonewtonian_coupled_var.Yield_Stress) 
    node.SetSolutionStepValue(INTERNAL_FRICTION_ANGLE,0,pfem_nonewtonian_coupled_var.Friction_Angle_Tangent)
##    node.SetSolutionStepValue(EFFECTIVE_VISCOSITY,0,pfem_nonewtonian_coupled_var.Yield_Stress)
##    node.SetSolutionStepValue(AUX_INDEX,0,0.0)

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
box_corner1[0]=pfem_nonewtonian_coupled_var.min_x;
box_corner1[1]=pfem_nonewtonian_coupled_var.min_y;
box_corner1[2]=pfem_nonewtonian_coupled_var.min_z;
box_corner2 = Vector(3);
box_corner2[0]=pfem_nonewtonian_coupled_var.max_x;
box_corner2[1]=pfem_nonewtonian_coupled_var.max_y;
box_corner2[2]=pfem_nonewtonian_coupled_var.max_z;



#the buffer size should be set up here after the mesh is read for the first time
pfem_model_part.SetBufferSize(2)
fixed_model_part.SetBufferSize(2)




fluid_viscosity   = edgebased_levelset_var.viscosity
##print "*****************   VISCOSITY  *********************"
fluid_density = edgebased_levelset_var.density
solid_density = pfem_nonewtonian_coupled_var.Density
##print "*****************   DENSITY  *********************"
fluid_solver = edgebased_levelset_solver.EdgeBasedLevelSetSolver(fixed_model_part,domain_size,body_force,fluid_viscosity,fluid_density)

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
#adding dofs to the structural solver
monolithic_solver_lagrangian_nonnewtonian.AddDofs(pfem_model_part)
#defining the structural solver
structural_solver = monolithic_solver_lagrangian_nonnewtonian.MonolithicSolver(pfem_model_part,domain_size,box_corner1,box_corner2)
#setting solver parameters
structural_solver.oss_swith = int(pfem_nonewtonian_coupled_var.use_oss)
structural_solver.dynamic_tau = pfem_nonewtonian_coupled_var.dynamic_tau
structural_solver.regularization_coef = pfem_nonewtonian_coupled_var.m_coef #m regularization coefficient in the exponential law of viscosity
structural_solver.echo_level = 2
#problem parameters to be changed
#####time setting pfem
### pfem implicit
max_dt = pfem_nonewtonian_coupled_var.max_dt
min_dt = pfem_nonewtonian_coupled_var.min_dt
### edgebased explicit
max_Dt = edgebased_levelset_var.max_time_step
initial_Dt_ls = 0.001 * max_Dt 

final_time = edgebased_levelset_var.max_time
output_dt = edgebased_levelset_var.output_dt
safety_factor = edgebased_levelset_var.safety_factor

number_of_inital_steps = edgebased_levelset_var.number_of_inital_steps
initial_time_step = edgebased_levelset_var.initial_time_step


structural_solver.Initialize(output_dt)



if (domain_size == 3):
    ProjectionModule = MeshTransfer3D();
else:
    ProjectionModule = MeshTransfer2D();
PfemUtils = PfemUtils();
DragUtils = DragUtils();
VariableUtils = VariableUtils()


ProjectionModule.DirectScalarVarInterpolation(pfem_model_part, fixed_model_part, POROSITY, POROSITY);
ProjectionModule.DirectScalarVarInterpolation(pfem_model_part, fixed_model_part, DIAMETER, DIAMETER);
for node in fixed_model_part.Nodes:
    if(node.GetSolutionStepValue(POROSITY) == 0.0):
        node.SetSolutionStepValue(POROSITY,0,1.0)
    if(node.GetSolutionStepValue(DIAMETER) == 0.0):
        node.SetSolutionStepValue(DIAMETER,0,1.0)



out = 0

original_max_dt = max_Dt
Time = 0.0
step = 0
next_output_time = output_dt
max_safety_factor = safety_factor

##substep_number = 1 #each substep_number structure is evaluated. TO BE PUT IN THE PROBLEM TYPE!!!!
##substep_counter = 0

check_structure_dt = True
print_structure = True
next_structure_step = Time

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

    ### CALCULATE FLUID_DT AND STRUCTURE_DT
    fluid_dt = fluid_solver.EstimateTimeStep(safety_factor,max_Dt)
    if(check_structure_dt == True):
        if(step < 2.0*number_of_inital_steps):
            structure_dt = min_dt
        else:
            structure_dt = structural_solver.EstimateDeltaTime(min_dt,max_dt)
        
        ##check which of the dime step is more restrictive
        if(structure_dt <= fluid_dt):
            Dt = structure_dt
            next_structure_step = Time
        else:
            Dt = fluid_dt
            if(structure_dt * safety_factor > 10.0* fluid_dt):
                next_structure_step = Time + 10.0* fluid_dt
            else:
                next_structure_step = Time + safety_factor * structure_dt
            check_structure_dt = False
    else:
        Dt = fluid_dt
                
    print "fluid_dt = ", fluid_dt
    print "structure_dt = ", structure_dt
    print "Dt = ", Dt
   
    print "next_structure_step =", next_structure_step

    Time = Time + Dt
    fixed_model_part.CloneTimeStep(Time)

    ##let's do this later ... only each substep_number fluid solutions
    ##pfem_model_part.CloneTimeStep(Time)

    print "******** CURRENT TIME = ",Time


            

    if(step > 3):
       
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


        if(Time >= next_structure_step):
            
            pfem_model_part.CloneTimeStep(Time)
            if(Time > 0.0):  ## only for A2!!!!
                
                       
                print pfem_model_part
                #solving the structural problem ...........................................................................................
                print "starting solving structure~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                PfemUtils.CalculateNodalMass(pfem_model_part,domain_size)

            

######                ##Projecting the water pressure
                ProjectionModule.DirectScalarVarInterpolation(fixed_model_part, pfem_model_part, PRESSURE, WATER_PRESSURE);


                #Calculating Darcy non linear term
                DragUtils.CalculateFluidDrag(fixed_model_part, SEEPAGE_DRAG, fluid_viscosity, fluid_density, solid_density )
                ProjectionModule.DirectVectorialVarInterpolation(fixed_model_part, pfem_model_part, SEEPAGE_DRAG, SEEPAGE_DRAG);
                ProjectionModule.DirectVectorialVarInterpolation(fixed_model_part, pfem_model_part, PRESS_PROJ, PRESS_PROJ);
                #Adding Darcy non linear term and Water weight to the external force                
                DragUtils.AddDrag(pfem_model_part, SEEPAGE_DRAG, BODY_FORCE, body_force )
                
                #PFEM steps
                structural_solver.Remesh();
                (structural_solver.solver).Solve();
                structural_solver.RestoreOldPosition()

                (structural_solver.PfemUtils).MoveLonelyNodes(pfem_model_part)
                (structural_solver.solver).Clear();
                
                check_structure_dt = True
                print_structure = True
                print "finished solving structure"
##        else:
##            for node in fixed_model_part.Nodes:
##                node.SetSolutionStepValue(DISPLACEMENT_X,0,0.0)
##                node.SetSolutionStepValue(DISPLACEMENT_Y,0,0.0)
##                node.SetSolutionStepValue(DISPLACEMENT_Z,0,0.0)
        

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
   	gid_io.WriteNodalResults(SEEPAGE_DRAG,fixed_model_part.Nodes,Time,0)


        gid_io.Flush()
        gid_io.FinalizeResults()

##        if(Time > 0.0):  ## only for ??? !!!!
        if(print_structure == True):  ## only for ??? !!!!
        
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
    ##        gid_io.WriteNodalResults(IS_BOUNDARY,pfem_model_part.Nodes,Time,0)
            gid_io.WriteNodalResults(IS_STRUCTURE,pfem_model_part.Nodes,Time,0)
            gid_io.WriteNodalResults(VISCOSITY,pfem_model_part.Nodes,Time,0)
            gid_io.WriteNodalResults(DENSITY,pfem_model_part.Nodes,Time,0)
            gid_io.WriteNodalResults(BODY_FORCE,pfem_model_part.Nodes,Time,0)
            gid_io.WriteNodalResults(DISPLACEMENT,pfem_model_part.Nodes,Time,0)
    ##        gid_io.WriteNodalResults(FRICTION_COEFFICIENT,pfem_model_part.Nodes,Time,0)
    ##        gid_io.WriteNodalResults(IS_EROSIONABLE,pfem_model_part.Nodes,Time,0)
            gid_io.WriteNodalResults(POROSITY,pfem_model_part.Nodes,Time,0)
            gid_io.WriteNodalResults(DIAMETER,pfem_model_part.Nodes,Time,0)
            gid_io.WriteNodalResults(SEEPAGE_DRAG,pfem_model_part.Nodes,Time,0)
            gid_io.WriteNodalResults(WATER_PRESSURE,pfem_model_part.Nodes,Time,0)
            gid_io.PrintOnGaussPoints(EQ_STRAIN_RATE,pfem_model_part,Time)
            gid_io.PrintOnGaussPoints(MU,pfem_model_part,Time)
            gid_io.PrintOnGaussPoints(TAU,pfem_model_part,Time)
            gid_io.WriteNodalResults(NODAL_H, pfem_model_part.Nodes, Time, 0);
            gid_io.WriteNodalResults(EXTERNAL_PRESSURE, pfem_model_part.Nodes, Time, 0);
            gid_io.WriteNodalResults(PRESS_PROJ, pfem_model_part.Nodes, Time, 0);
            

            gid_io.Flush()  
            gid_io.FinalizeResults()

            print_structure = False
            
        next_output_time = Time + output_dt
        
        out = 0

    out = out + 1
    step = step + 1
    print "step finished"

print "solution finished"        
