 # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
 # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import fluid_only_var

#
#
print("Click2cast running state: Start KRATOS")


#
#
# ATTENTION: here the order is important

# from now on the order is not anymore crucial
#
import sys
import math
# import time as timer
import time as time_measure
# sys.path.append(fluid_only_var.kratos_path)
sys.argv.append("ProjectParameters.conf")
print(sys.argv)
if len(sys.argv) < 2:
    print ("ERROR: The configuration file is missing")
    print ("")
    print ("Syntax:")
    print ("    KratosFastFill.exe configfile.conf")
    print ("")
    print ("example:")
    print ("    KratosFastFill.exe ProjectParameters.conf")
    sys.exit()
	
# including kratos path
exec(open(sys.argv[1]).read())

# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

start_time = time_measure.time()

# importing Kratos main library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.Click2CastApplication import *

verifyp = Click2CastVerifyProcess(sys.argv[0])
verifyp.Execute()


# assign inlet velocity

# def ComputeAreaToNode(reference_node,ThisModelPart):
    # area=Vector(3)
    # for condition in ThisModelPart.Conditions:
        # for node in condition.GetNodes():
            # if node.Id==reference_node.Id:
                # normal = condition.GetValue(NORMAL)
                # #normal_size = math.sqrt(normal * normal)
                # area=area + normal #_size
    # return area
        
# def ApplyInlet(InletConditions, ThisModelPart):
    # inlet_velocity = 0.00
    # current_time = 0.00
    # #Set velocities to 0.
    # for condition in InletConditions:
        
        # for node in condition.GetNodes():
            # velocity = Vector(3)
            # velocity[0]=0.0
            # velocity[1]=0.0
            # velocity[2]=0.0
            # node.SetSolutionStepValue(VELOCITY, velocity)
            
    # if(hasattr(ProjectParameters, "InletTableId")):
        # inlet_velocity = ThisModelPart.GetTable(
            # ProjectParameters.InletTableId).GetValue(current_time)
    # else:
        # inlet_velocity = ProjectParameters.inlet_velocity
    # for condition in InletConditions:
        # normal = condition.GetValue(NORMAL)
        # normal_size = math.sqrt(normal * normal)

        # for node in condition.GetNodes():
            # print("node ____")
            # print(node.Id)
            # inc_velocity=node.GetSolutionStepValue(VELOCITY)
            # #print(ComputeAreaToNode(node,ThisModelPart))
            # projection=1.0/(ComputeAreaToNode(node,ThisModelPart)*normal/normal_size)
            # velocity = inc_velocity+(inlet_velocity *  (normal * projection))
            # node.SetSolutionStepValue(VELOCITY, velocity)
            # node.SetSolutionStepValue(
                # TEMPERATURE, ProjectParameters.FLUID_TEMPERATURE)
            # node.Fix(VELOCITY_X)
            # node.Fix(VELOCITY_Y)
            # node.Fix(VELOCITY_Z)
            # node.Fix(TEMPERATURE)
            # # node.SetSolutionStepValue(DISTANCE,0,-1.0)
            # # node.Fix(DISTANCE)
            # node.SetValue(Y_WALL, 0.000)
            # node.SetSolutionStepValue(IS_SLIP, 0, 0.0)
            # node.SetSolutionStepValue(IS_STRUCTURE, 0, 0.0)
            # node.Set(INLET)

            # node.SetValue(IS_STRUCTURE, 0.0)
    # return inlet_velocity

def ApplyInlet(InletConditions, ThisModelPart):
    inlet_velocity = 0.00
    current_time = 0.00
    if(hasattr(ProjectParameters, "InletTableId")):
        inlet_velocity = ThisModelPart.GetTable(
            ProjectParameters.InletTableId).GetValue(current_time)
    else:
        inlet_velocity = ProjectParameters.inlet_velocity
    for condition in InletConditions:
        normal = condition.GetValue(NORMAL)
        normal_size = math.sqrt(normal * normal)

        for node in condition.GetNodes():
            velocity = 1.0 * inlet_velocity / normal_size * \
                normal
            node.SetSolutionStepValue(VELOCITY, velocity)
            node.SetSolutionStepValue(
                TEMPERATURE, ProjectParameters.FLUID_TEMPERATURE)
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.Fix(VELOCITY_Z)
            node.Fix(TEMPERATURE)
            node.SetValue(Y_WALL, 0.000)
            node.SetSolutionStepValue(IS_SLIP, 0, 0.0)
            node.SetSolutionStepValue(IS_STRUCTURE, 0, 0.0)
            node.Set(INLET)

            node.SetValue(IS_STRUCTURE, 0.0)
    return inlet_velocity
    


def ApplyInletVelocity(InletConditions, ThisModelPart, Step, MaxEdge):
    inlet_velocity = 0.00
    too_low = 0
    zero_vel = Vector(3)
    zero_vel[0] = 0.0
    zero_vel[1] = 0.0
    zero_vel[2] = 0.0
    current_time = ThisModelPart.ProcessInfo[TIME]
    if(hasattr(ProjectParameters, "InletTableId")):
        inlet_pressure = ThisModelPart.GetTable(
            ProjectParameters.InletTableId).GetValue(current_time)
        print (inlet_pressure)
        inlet_velocity = math.sqrt(2.0 * inlet_pressure)
        inlet_node = InletConditions[0].GetNodes()[0]
# print(inlet_node)
        print(inlet_velocity)
        surfaceNode = Vector(3)
        maxh = -10000000.0
        for nd in ThisModelPart.Nodes:
            dist = nd.GetSolutionStepValue(DISTANCE)
            if (math.fabs(dist) < MaxEdge):
                xx = nd.X - inlet_node.X
                yy = nd.Y - inlet_node.Y
                zz = nd.Z - inlet_node.Z
                distoinlet = xx * xx + yy * yy + zz * zz
                if(distoinlet > maxh):
                    maxh = distoinlet
                    surfaceNode[0] = xx
                    surfaceNode[1] = yy
                    surfaceNode[2] = zz

        Effective_height = -(surfaceNode[0] * ProjectParameters.body_force_x + surfaceNode[1] * ProjectParameters.body_force_y +
                             surfaceNode[2] * ProjectParameters.body_force_z)
        print("''''''''''''''''''''''")
        print(Effective_height)
        print("''''''''''''''''''''''")
        if(Effective_height > 0.0):
            inlet_pressure -= Effective_height
            if(inlet_pressure <= 0.0):
                print (
                    ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>******** 00000000000000000000000000000 *************")
                print (
                    ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>********* STOP DUE TO LOW INPUT PRESSURE ***************")
                print (
                    ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>******** 00000000000000000000000000000 *************")
                return 1.0
            inlet_velocity = math.sqrt(2.0 * inlet_pressure)
            print(inlet_velocity)
    else:
        inlet_velocity = ProjectParameters.inlet_velocity
    if(hasattr(ProjectParameters, "second_phase_node")):
        if (ThisModelPart.Nodes[ProjectParameters.second_phase_node].GetSolutionStepValue(DISTANCE) < 0.00):
            inlet_velocity = ProjectParameters.second_phase_inlet_velocity
            # for the next iterations to be always used
            ProjectParameters.inlet_velocity = inlet_velocity

    for condition in InletConditions:
        normal = condition.GetValue(NORMAL)
        normal_size = math.sqrt(normal * normal)
        for node in condition.GetNodes():
            velocity = 1.0 * inlet_velocity / normal_size * \
                normal
            node.SetSolutionStepValue(VELOCITY, velocity)
            if(step < 3 or 1 == 1):
                node.SetSolutionStepValue(VELOCITY, 1, velocity)
                node.SetSolutionStepValue(VELOCITY, 2, velocity)
                node.SetSolutionStepValue(MESH_VELOCITY, 1, zero_vel)
                node.SetSolutionStepValue(MESH_VELOCITY, 2, zero_vel)
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.Fix(VELOCITY_Z)
    return too_low


def IncreaseWallLawInSolidifiedZone(ThisModelPart, Ywall):
    for node in ThisModelPart.Nodes:
        slip_flag = node.GetSolutionStepValue(IS_SLIP)
        distance = node.GetSolutionStepValue(DISTANCE)
        if(slip_flag != 0.0 and distance<0.0):
            alpha = node.GetSolutionStepValue(SOLID_FRACTION)
            
            modification_coefficient = 1.0+100.0*alpha*alpha
            if(alpha != 0.0):
                node.SetValue(Y_WALL, modification_coefficient * Ywall)

timer = Timer()

# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")
thermal_model_part = ModelPart("ThermalPart")
solidification_model_part = ModelPart("SolidificationPart")
#
# importing the solvers needed
import dpg_monolithic_solver_eulerian
dpg_monolithic_solver_eulerian.AddVariables(fluid_model_part)

# For thermal analysis
my_settings = ConvectionDiffusionSettings()

my_settings.SetDensityVariable(DENSITY)
my_settings.SetDiffusionVariable(CONDUCTIVITY)
my_settings.SetUnknownVariable(TEMPERATURE)
my_settings.SetVolumeSourceVariable(HEAT_FLUX)
my_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
my_settings.SetMeshVelocityVariable(MESH_VELOCITY)
my_settings.SetConvectionVariable(VELOCITY)
my_settings.SetTransferCoefficientVariable(HTC)
import convdiff_phasechange_solver

convdiff_phasechange_solver.AddVariables(fluid_model_part, my_settings)
convdiff_phasechange_solver.AddVariables(thermal_model_part, my_settings)
convdiff_phasechange_solver.AddVariables(solidification_model_part, my_settings)
fluid_model_part.AddNodalSolutionStepVariable(CONDUCTIVITY)
fluid_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
fluid_model_part.AddNodalSolutionStepVariable(HEAT_FLUX)
fluid_model_part.AddNodalSolutionStepVariable(FACE_HEAT_FLUX)
fluid_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
fluid_model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT)
fluid_model_part.AddNodalSolutionStepVariable(TEMPERATURES)
fluid_model_part.AddNodalSolutionStepVariable(NODAL_MASS)
fluid_model_part.AddNodalSolutionStepVariable(FRONT_MEETING)
fluid_model_part.AddNodalSolutionStepVariable(MACRO_POROSITY)
fluid_model_part.AddNodalSolutionStepVariable(SHRINKAGE_POROSITY)
fluid_model_part.AddNodalSolutionStepVariable(FILLTIME)
fluid_model_part.AddNodalSolutionStepVariable(REACTION)
if(hasattr(ProjectParameters, "SI_UNITS")):
    if(ProjectParameters.SI_UNITS == 0):
        fluid_model_part.AddNodalSolutionStepVariable(SHRINKAGE_POROSITY_US)
        fluid_model_part.AddNodalSolutionStepVariable(SOLIDIF_MODULUS_US)
        fluid_model_part.AddNodalSolutionStepVariable(TEMPERATURES_US)

# introducing input file name
input_file_name = ProjectParameters.problem_name
log_file = open(input_file_name + ".log", 'w')
time_log_file = open(input_file_name + ".tlog", 'w')
solver_log_file = open(input_file_name + ".slog", 'w')


timer.Start("Reading Model")

# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly  # WriteConditions
# write_conditions)
gid_io = GidIO(input_file_name + "_F_k", gid_mode,
               multifile, deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)
timer.Stop("Reading Model")
timer.Start("Initialization")

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)
thermal_model_part.SetBufferSize(3)
# adding dofs
dpg_monolithic_solver_eulerian.AddDofs(fluid_model_part)
convdiff_phasechange_solver.AddDofs(fluid_model_part, my_settings)

# copy to vulcan!!!!!!!!!!!!!!!!!
copy_filling_variables_process = DPGCopyToVulcanPostVariablesProcess(
    fluid_model_part, 1, 0, 1)
copy_solidification_variables_process = DPGCopyToVulcanPostVariablesProcess(
    fluid_model_part, 0, 1, 1)
if(hasattr(ProjectParameters, "SI_UNITS")):
    if(ProjectParameters.SI_UNITS == 1):
        copy_filling_variables_process = DPGCopyToVulcanPostVariablesProcess(
            fluid_model_part, 1, 0, 1)
        copy_solidification_variables_process = DPGCopyToVulcanPostVariablesProcess(
            fluid_model_part, 0, 1, 1)
    else:
        copy_filling_variables_process = DPGCopyToVulcanPostVariablesProcess(
            fluid_model_part, 1, 0, 0)
        copy_solidification_variables_process = DPGCopyToVulcanPostVariablesProcess(
            fluid_model_part, 0, 1, 0)

# creating the solvers
# fluid solver wall_law_y
mat_density = (fluid_model_part.GetTable(1)).GetValue(
    ProjectParameters.FLUID_TEMPERATURE)
if(ProjectParameters.max_time > 1.0):  # gravity_filling
    y_wall_val = 10.0 * ProjectParameters.wall_law_y
    y_wall_fac = 1.0  # 1.0
    C_SMAG = 0.45  # 45#ProjectParameters.SmagorinskyConstant
    CFL = 10.0
    max_itr_solver = 4  # 8
    volume_correction_switch = True
    air_exit_flag = False
    is_gravity_filling = 0
    corrected_rho2 = 1.0  # 0.04*mat_density
    mold_temp = ProjectParameters.AMBIENT_TEMPERATURE
else:
    y_wall_val = ProjectParameters.wall_law_y
    y_wall_fac = 10000.0
    C_SMAG = 0.25  # ProjectParameters.SmagorinskyConstant
    CFL = 10.0
    max_itr_solver = 4
    volume_correction_switch = True
    air_exit_flag = False
    is_gravity_filling = 0
    corrected_rho2 = 1.0
    mold_temp = ProjectParameters.AMBIENT_TEMPERATURE

max_itr_solver = 1#it does not really improve with more iterations!!

# writing log
log_file.write("y_wall_val : " + str(y_wall_val) + "\n")
log_file.write("y_wall_fac : " + str(y_wall_fac) + "\n")
log_file.write("C_SMAG     : " + str(C_SMAG) + "\n")
log_file.write("CFL        : " + str(CFL) + "\n")
log_file.write("max_itr_solver           : " + str(max_itr_solver) + "\n")
log_file.write("volume_correction_switch : " +
               str(volume_correction_switch) + "\n")
log_file.write("air_exit_flag            : " + str(air_exit_flag) + "\n")
log_file.write("is_gravity_filling       : " + str(is_gravity_filling) + "\n")
log_file.write("corrected_rho2           : " + str(corrected_rho2) + "\n")


# Filling and solidification switches
FILLING = ProjectParameters.filling
SOLIDIFICATION = ProjectParameters.solidification

air_temp = ProjectParameters.FLUID_TEMPERATURE 
initial_density = fluid_model_part.GetTable(1).GetValue(ProjectParameters.AMBIENT_TEMPERATURE)
for node in fluid_model_part.Nodes:
    node.SetValue(Y_WALL, y_wall_val)
    node.SetSolutionStepValue(DISTANCE, 0, 1000.0)
    node.SetSolutionStepValue(DENSITY, 0, initial_density)
    node.SetSolutionStepValue(VELOCITY_X, 0, 0.0)
    node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)
    node.SetSolutionStepValue(VELOCITY_Z, 0, 0.0)
    node.SetSolutionStepValue(BODY_FORCE_X, 0, ProjectParameters.body_force_x)
    node.SetSolutionStepValue(BODY_FORCE_Y, 0, ProjectParameters.body_force_y)
    node.SetSolutionStepValue(BODY_FORCE_Z, 0, ProjectParameters.body_force_z)
    node.SetSolutionStepValue(
        TEMPERATURE, 0, air_temp) #ProjectParameters.SOLID_TEMPERATURE)
    node.SetSolutionStepValue(
        TEMPERATURE, 1, air_temp) #ProjectParameters.SOLID_TEMPERATURE)
    node.SetSolutionStepValue(
        TEMPERATURE, 2, air_temp) #ProjectParameters.SOLID_TEMPERATURE)
for condition in fluid_model_part.Conditions:
    # To compute Normals and slip
    for node in condition.GetNodes():
        node.SetSolutionStepValue(IS_SLIP, 0, 1.0)
        node.SetSolutionStepValue(IS_STRUCTURE, 0, 1.0)
        node.SetValue(IS_STRUCTURE, 1.0)
        if(condition.GetValue(IS_INLET) > 0.0):
            node.SetSolutionStepValue(DISTANCE, 0, -100.0)
            node.Fix(DISTANCE)

fluid_solver = dpg_monolithic_solver_eulerian.MonolithicSolver(
    fluid_model_part, domain_size)
oss_swith = 0
dynamic_tau = 0.0 
fluid_solver.oss_switch = oss_swith
fluid_solver.dynamic_tau_fluid = dynamic_tau
fluid_solver.dynamic_tau_levelset = 0.001
fluid_solver.max_iter = max_itr_solver
fluid_solver.echo_level = 0
fluid_solver.CalculateReactionFlag = False
fluid_solver.ReformDofSetAtEachStep = False
fluid_solver.use_slip_conditions = True
fluid_solver.volume_correction_switch = volume_correction_switch
fluid_solver.redistance_frequency = 1
fluid_solver.CFL = CFL
fluid_solver.vol_cr_step = 0
fluid_solver.rho1 = mat_density
fluid_solver.rho2 = corrected_rho2
fluid_solver.rel_vel_tol = 3e-2
fluid_solver.rel_pres_tol = 3e-2
fluid_solver.Initialize()

# generate temperature model part
modeler_util = ConnectivityPreserveModeler()
print (thermal_model_part)

if(domain_size == 2):
    modeler_util.GenerateModelPart(
        fluid_model_part, thermal_model_part, "SUPGConvDiff2D", "Condition2D")
else:
    modeler_util.GenerateModelPart(fluid_model_part, thermal_model_part, "SUPGConvDiffPhaseChange3DLinearized", "VirtualMouldElement3D")
    for cond in fluid_model_part.Conditions:
        if (cond.GetValue(IS_STRUCTURE) == 1):
            for node in cond.GetNodes():
                node.SetSolutionStepValue(IS_BOUNDARY, 0, 1.0)
				

fluid_model_part.ProcessInfo.SetValue(
    AMBIENT_TEMPERATURE, mold_temp)
fluid_model_part.ProcessInfo.SetValue(
    SOLID_TEMPERATURE, ProjectParameters.SOLID_TEMPERATURE)
fluid_model_part.ProcessInfo.SetValue(
    LATENT_HEAT, ProjectParameters.LATENT_HEAT)

fluid_model_part.ProcessInfo.SetValue(
    DENSITY, (fluid_model_part.GetTable(1)).GetValue(ProjectParameters.AMBIENT_TEMPERATURE))
fluid_model_part.ProcessInfo.SetValue(
    SPECIFIC_HEAT, (fluid_model_part.GetTable(2)).GetValue(ProjectParameters.AMBIENT_TEMPERATURE))
fluid_model_part.ProcessInfo.SetValue(
    HTC, (fluid_model_part.GetTable(7)).GetValue(ProjectParameters.AMBIENT_TEMPERATURE))
fluid_model_part.ProcessInfo.SetValue(K0, ProjectParameters.MODULUS)
fluid_model_part.ProcessInfo.SetValue(IS_GRAVITY_FILLING, is_gravity_filling)

# Now we add the variables needed for the Virtual Mould
fluid_model_part.ProcessInfo.SetValue( MOULD_DENSITY, ProjectParameters.MOULD_DENSITY)
fluid_model_part.ProcessInfo.SetValue( MOULD_SPECIFIC_HEAT, ProjectParameters.MOULD_SPECIFIC_HEAT)
fluid_model_part.ProcessInfo.SetValue( MOULD_THICKNESS, ProjectParameters.MOULD_THICKNESS)
fluid_model_part.ProcessInfo.SetValue( MOULD_SFACT, ProjectParameters.MOULD_SFACT)
fluid_model_part.ProcessInfo.SetValue( MOULD_VFACT, ProjectParameters.MOULD_VFACT)
fluid_model_part.ProcessInfo.SetValue( MOULD_CONDUCTIVITY, ProjectParameters.MOULD_CONDUCTIVITY)
fluid_model_part.ProcessInfo.SetValue( MOULD_HTC_ENVIRONMENT, ProjectParameters.MOULD_HTC_ENVIRONMENT)

them_solver = convdiff_phasechange_solver.Solver(
    thermal_model_part, domain_size, my_settings)
them_solver.max_distance = 0.9*fluid_solver.max_distance
them_solver.MaxLineSearchIterations = 3
them_solver.Initialize()

inlet_conditions = []
for condition in fluid_model_part.Conditions:
    if(condition.GetValue(IS_INLET) > 0.0):
        inlet_conditions.append(condition)

inlet_vel = ApplyInlet(inlet_conditions, fluid_model_part)

input_discharge = ProjectParameters.input_discharge
fill_time = ProjectParameters.max_time
tot_volume = input_discharge * fill_time
max_Dt = 5.0 * ProjectParameters.max_time_step
nominal_inlet_vel = ProjectParameters.inlet_velocity
inlet_area = input_discharge / nominal_inlet_vel
input_discharge = inlet_vel * inlet_area
print("Fill Time")
print(fill_time)
#fill_time = tot_volume / input_discharge
max_Dt = fill_time / 200.0



full_Dt = max_Dt
initial_Dt = 1.0 * full_Dt  # 0.05 #0.01
output_step = 1.0

out = output_step
output_dt = ProjectParameters.output_dt
next_output_time = output_dt
next_output_time = 1e-6


##inlet_vel = ProjectParameters.inlet_velocity

# Smagorinsky
print ("The C_SMAG is: ", C_SMAG)
for elem in fluid_model_part.Elements:
    elem.SetValue(C_SMAGORINSKY, C_SMAG)

time = 0.0
step = 0

def OpenAirExitInFarDryZone(model_part, bx, by, bz):
    body_force = Vector(3)
    body_force[0] = bx
    body_force[1] = by
    body_force[2] = bz
    zero = Vector(3)
    zero[0] = 0.0
    zero[1] = 0.0
    zero[2] = 0.0

    number_of_exits = 0
    max_d = 0.0
    for node in model_part.Nodes:
        d = node.GetSolutionStepValue(DISTANCE)
        if(d > max_d):
            max_d = d

    out_skin_list = []

    for node in model_part.Nodes:
        slip_flag = node.GetSolutionStepValue(IS_SLIP)
        dist = node.GetSolutionStepValue(DISTANCE)
        node.SetSolutionStepValue(IS_FREE_SURFACE, 0, 0.0)
        if(dist >= 0.9 * max_d):
            node.SetSolutionStepValue(BODY_FORCE, 0, zero)
            node.SetSolutionStepValue(IS_FREE_SURFACE, 0, 10.0)
            node.SetSolutionStepValue(VELOCITY,0,zero)
            node.SetSolutionStepValue(VELOCITY,1,zero)
            node.SetSolutionStepValue(VELOCITY,2,zero)
            number_of_exits += 1

            node.Fix(PRESSURE)
            node.SetSolutionStepValue(PRESSURE, 0, 0.0)
            if (slip_flag != 0):

                node.SetValue(IS_STRUCTURE, 0.0)
        else:
            node.SetSolutionStepValue(BODY_FORCE, 0, body_force)
            node.Free(PRESSURE)
            if (slip_flag == 10.0 or slip_flag == 20.0 or slip_flag == 30.0):

                if(dist > 0.0):
                    if(slip_flag != 10):
                        
                        vel = node.GetSolutionStepValue(VELOCITY)
                        n = node.GetSolutionStepValue(NORMAL)
                        prod = n[0]*vel[0] + n[1]*vel[1] + n[1]*vel[1]
                        if(prod > 0.0): ##going outside!
                            
                            node.Fix(PRESSURE)
                            node.SetSolutionStepValue(VELOCITY,0,zero)   #newly added
                            node.SetSolutionStepValue(VELOCITY,1,zero)   #newly added
                            node.SetSolutionStepValue(VELOCITY,2,zero)   #newly added
                            node.SetSolutionStepValue(PRESSURE, 0, 0.0)
                            node.SetValue(IS_STRUCTURE, 0.0)
                            node.SetValue(Y_WALL, 100.0 * y_wall_val)
                        else:                             
                            node.SetValue(IS_STRUCTURE, 1.0)
                    else:
                        out_skin_list.append(node)
                        node.SetValue(IS_STRUCTURE, 1.0)
                if(dist < 0.0):  # and slip_flag == 30.0
                    node.SetValue(IS_STRUCTURE, 1.0)
                    if(slip_flag != 10):
                        node.SetValue(Y_WALL, 100.0 * y_wall_val)

    if(number_of_exits == 0 and len(out_skin_list) != 0):
        max_d = 0.0
        for node in out_skin_list:
            d = node.GetSolutionStepValue(DISTANCE)
            if(d > max_d):
                max_d = d

        dlimit = 0.5 * max_d

        for node in out_skin_list:
            d = node.GetSolutionStepValue(DISTANCE)
            if(d > dlimit):
                node.Fix(PRESSURE)
                node.SetSolutionStepValue(PRESSURE, 0, 0.0)
                node.SetValue(IS_STRUCTURE, 0.0)
                node.SetSolutionStepValue(IS_FREE_SURFACE, 0, 20.0)

all_exits_blocked = 1
temp_time = 0.0
step = 0
last_print = -1
is_dry = 1.0
is_hot = 1.0
too_low = 0.0
screen_output_dt = 1
next_screen_output = screen_output_dt
time_output_dt = 10.00 * screen_output_dt
next_time_output = time_output_dt
filled_percent = 0.0
net_volume = 0.0
fluid_model_part.ProcessInfo.SetValue(NET_INPUT_MATERIAL, 0.0)
write_US_units = False
if(hasattr(ProjectParameters, "SI_UNITS")):
    if(ProjectParameters.SI_UNITS == 0):
        write_US_units = True


timer.Stop("Initialization")

import thermal_solution_utilities
thermal_utilities = thermal_solution_utilities.ThermalSolutionUtilities(thermal_model_part,ProjectParameters)

if(FILLING == 1.0):
    print("Click2cast running state: Start FILLING solver")
    sys.stdout.flush()
    # mesh to be printed
    mesh_name = 0.0
    gid_io.InitializeMesh(mesh_name)
    gid_io.WriteMesh(fluid_model_part.GetMesh())
    gid_io.FinalizeMesh()
    gid_io.Flush()
    gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())

    air_entrapment_recognition_process = AirEntrapmentRecognitionProcess(
        fluid_model_part)
    front_meeting_recognition_process = FrontMeetingRecognitionProcess(
        fluid_model_part)
    if(hasattr(ProjectParameters, "TiltPouring")):
        body_force = Vector(3)
        body_force[0] = ProjectParameters.body_force_x
        body_force[1] = ProjectParameters.body_force_y
        body_force[2] = ProjectParameters.body_force_z
        center_point = Vector(3)
        center_point[0] = ProjectParameters.RotationPointX
        center_point[1] = ProjectParameters.RotationPointY
        center_point[2] = ProjectParameters.RotationPointZ

        tilt_pouring_process = TiltPouringProcess(
            fluid_model_part, ProjectParameters.InitialAngle, ProjectParameters.max_time, center_point, ProjectParameters.RotationAxis, body_force)
    test = open("test.csv", 'w')

    time = time + max_Dt * 0.0001
    fluid_model_part.CloneTimeStep(time)
    
    done_once = False

    while(filled_percent <= 99.5 and is_dry == 1.0 and is_hot == 1.0 and too_low == 0.0):

        if(hasattr(ProjectParameters, "TiltPouring")):
            tilt_pouring_process.Execute()
        #percent_done = 100.00 * (time / final_time)
        Dt = EstimateTimeStep3D().ComputeDt(
            fluid_model_part, 4.0 * fluid_solver.max_edge_size, CFL, 0.05 * max_Dt, max_Dt)
        time = time + Dt
        time_crt = time
        fluid_model_part.CloneTimeStep(time)
        # elapsed_tiirme = timer.clock() - start_time
        elapsed_time = time_measure.time() - start_time
        print ("DT= ", Dt, "max_dt= ", max_Dt, "Total time= ", elapsed_time)

        # thermal copy
        thermal_model_part.ProcessInfo = fluid_model_part.ProcessInfo
         # thermal solve
        thermal_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.0)
        thermal_utilities.BeforeThermalSolution(fluid_solver.max_distance)


        #thermal_model_part.ProcessInfo.SetValue(LATENT_HEAT,0.0)
        them_solver.Solve()  
        
        
        BiphasicFillingUtilities().ApplyTemperatureLimitation(fluid_model_part,ProjectParameters.FLUID_TEMPERATURE, ProjectParameters.AMBIENT_TEMPERATURE)
        
        thermal_utilities.AfterThermalSolution()

        is_hot = BiphasicFillingUtilities().SolidificationDuringFilling(fluid_model_part, fluid_solver.max_edge_size)

        if (is_hot == 0):
            print (
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>******** 00000000000000000000000000000 *************")
            print (
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>********* STOP DUE TO SOLIDIFICATION ***************")
            print (
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>******** 00000000000000000000000000000 *************")        
        

        BiphasicFillingUtilities().ComputeNetInletVolume(fluid_model_part)
        timer.Start("Fluid Solve")
        gravity = Vector(3)
        if(hasattr(ProjectParameters, "TiltPouring")):
            for node in fluid_model_part.Nodes:
                gravity = node.GetSolutionStepValue(BODY_FORCE)
                break
        else:
            gravity[0] = ProjectParameters.body_force_x
            gravity[1] = ProjectParameters.body_force_y
            gravity[2] = ProjectParameters.body_force_z
         
        if(done_once == True):
            ActivationUtilities().ActivateElementsAndConditions( fluid_model_part, DISTANCE, 0.9*fluid_solver.max_distance, True) 
        done_once = True
        
        OpenAirExitInFarDryZone(
            fluid_model_part, gravity[0], gravity[1], gravity[2])
        too_low = ApplyInletVelocity(
            inlet_conditions, fluid_model_part, step, fluid_solver.maxmin[0])
        
        #Open air escape after all edge exits are closed
        
        fluid_solver.Solve(step)
        
        max_acceptable_acc_norm = 6.0 * ProjectParameters.inlet_velocity / max_Dt
        BiphasicFillingUtilities().ApplyVelocityLimitation(
            fluid_model_part, max_acceptable_acc_norm)

        solver_log_file.write("Step : " + str(step) + "\t: ")
        solver_log_file.write(str(fluid_solver.linear_solver) + "\n")

        timer.Stop("Fluid Solve")

        # correct time
        time_fact = 1.0
        if(volume_correction_switch == True and step > fluid_solver.vol_cr_step):
            wet_vol = fluid_model_part.ProcessInfo[WET_VOLUME]
            time_fact = wet_vol / tot_volume
            #time_crt = fill_time * time_fact
            time_crt = fill_time*time_fact
            print("//// Corrections: wet_vol = %e" % wet_vol, " tot_volume= %e" % tot_volume)
            print("//// Corrections: time = %e" % time, " time_crt= %e" % time_crt, " time_fact = %e" % time_fact)

        timer.Start("Air entrapment")
        # Air entrapment
        air_entrapment_recognition_process.Execute()
        timer.Stop("Air entrapment")
        timer.Start("Front meeting")
        # Front meeting
        front_meeting_recognition_process.Execute()
        timer.Stop("Front meeting")
        # Decide air escape
        # filled_percent=BiphasicFillingUtilities().AssignSmoothBoundaryAirExit(fluid_model_part,air_exit_flag,y_wall_val,y_wall_fac)
        filled_percent = BiphasicFillingUtilities().ComputeFillPercentage(
            fluid_model_part, time_crt)
 
        ## Open air escape after all edge exits are closed
                 
        if(int(filled_percent) >= int(next_screen_output)):
            print ()
            print ("Filled %.0f" % filled_percent, "% \t", "%e" %
                   time, "\t", "%e" % Dt, "\t", elapsed_time)
            log_file.write("Step : " + str(step) + "\t Filled {:.1f}".format(
                filled_percent) + "\t {:e}".format(time) + "\t {:e}".format(Dt) + "\t" + str(elapsed_time) + "\n")
            sys.stdout.flush()
            while(next_screen_output < filled_percent):
                next_screen_output += screen_output_dt
                print ("next output", next_screen_output)
            # print Timer()
        aaa = "step= " + str(step) + " time= " + str(
            time) + " time_crt= " + str(time_crt) + "\n"
        test.write(aaa)
        

                
        # if(out == output_step):
        if(time_crt >= next_output_time):
            verifyp.Execute()
            timer.Start("Write Results")
                
            timer.Start("Copy IO variables")
            copy_filling_variables_process.Execute()
            
            timer.Stop("Copy IO variables")
            if(hasattr(ProjectParameters, "TiltPouring")):
                tilt_pouring_process.ExecuteBeforeOutputStep()              

            gid_io.WriteNodalResults(
                VELOCITY, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                DISTANCE, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                MATERIAL, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                TEMPERATURE, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                TEMPERATURES, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                LAST_AIR, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                DENSITY, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                ENTHALPY, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                VELOCITIES, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                FRONT_MEETING, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                PRESSURES, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                FILLTIME, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                SOLID_FRACTION, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                DISPLACEMENT, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                MAX_VEL, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                HTC, fluid_model_part.Nodes, time_crt, 0)
            if(write_US_units):
                gid_io.WriteNodalResults(
                    TEMPERATURES_US, fluid_model_part.Nodes, time_crt, 0)
            gid_io.WriteNodalResults(
                PRESSURE, fluid_model_part.Nodes, time_crt, 0)

            last_print = time
            gid_io.Flush()
            sys.stdout.flush()

            out = 0
            next_output_time += output_dt
            timer.Stop("Write Results")
        if(int(filled_percent) >= next_time_output):
            print (timer)
            time_log_file.write("Step : " + str(step) + "\t Filled {:.1f}".format(
                filled_percent) + "\t {:e}".format(time) + "\t {:e}".format(Dt) + "\t" + str(elapsed_time) + "\n")
            time_log_file.write(str(timer) + "\n")
            next_time_output += time_output_dt
            if(next_time_output < time_crt + output_dt):
                next_time_output = time_crt + output_dt

        out = out + 1
        step = step + 1
    # write the last step results
    is_dry = BiphasicFillingUtilities().CheckIfAllNodesAreWet(fluid_model_part)
    if((last_print != time or is_dry == 1.0)and is_hot == 1.0 and too_low == 0.0):
        print("**************************LAST PRINT********************************")
        time_crt *= 1.01
        BiphasicFillingUtilities().LastStepExtrapolations(
            fluid_model_part, time_crt)
        copy_filling_variables_process.Execute()
        gid_io.WriteNodalResults(MATERIAL, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(DISTANCE, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(
            TEMPERATURES, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(
            VELOCITIES, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(FILLTIME, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(
            FRONT_MEETING, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(
                ENTHALPY, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(
            SOLID_FRACTION, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(
            DISPLACEMENT, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(LAST_AIR, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(MAX_VEL, fluid_model_part.Nodes, time_crt, 0)
        if(write_US_units):
            gid_io.WriteNodalResults(
                TEMPERATURES_US, fluid_model_part.Nodes, time_crt, 0)
        gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time_crt, 0)

        gid_io.FinalizeResults()
        gid_io.Flush()
        sys.stdout.flush()
    gid_io.FinalizeResults()
    print("Click2cast running state: End FILLING solver")
    # Write last temperature
    Last_temperature = open(input_file_name + ".Tmp", 'w')
    for node in fluid_model_part.Nodes:
        Last_temperature.write(str(node.Id) + "\t" + str(
            node.GetSolutionStepValue(TEMPERATURE)) + "\n")
    
    gid_io.Flush()
    sys.stdout.flush()
    Last_temperature.close()
    gid_io.Flush()
    sys.stdout.flush()

    
gid_io = 0


if(SOLIDIFICATION == 1.0 and FILLING==0):
    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(DENSITY, 0, initial_density)
        T = ProjectParameters.FLUID_TEMPERATURE
        node.SetSolutionStepValue(TEMPERATURE,0,T)
        node.SetSolutionStepValue(TEMPERATURE,1,T)
        node.SetSolutionStepValue(TEMPERATURE,2,T)
 

if(SOLIDIFICATION == 1.0):
    verifyp.Execute()
    print ("Click2cast running state: Start SOLIDIFICATION-COOLING solver")
    sys.stdout.flush()
    gid_io = GidIO(input_file_name + "_S_k", gid_mode,
                   multifile, deformed_mesh_flag, write_conditions)

    fluid_model_part.ProcessInfo.SetValue(LATENT_HEAT, ProjectParameters.LATENT_HEAT)

    # mesh to be printed
    mesh_name = 0.0
    gid_io.InitializeMesh(mesh_name)
    gid_io.WriteMesh(fluid_model_part.GetMesh())
    gid_io.FinalizeMesh()
    gid_io.Flush()
    gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())
    solidification_model_part.SetBufferSize(3)

    modeler_util.GenerateModelPart(fluid_model_part, solidification_model_part,
                                   "SUPGConvDiffPhaseChange3D", "EnvironmentContact3D")  # SUPGConvDiffPhaseChange3D
    for cond in fluid_model_part.Conditions:
        for node in cond.GetNodes():
            node.SetSolutionStepValue(IS_BOUNDARY, 0, 1.0)
            node.SetSolutionStepValue(DISTANCE, 0, -1.0)
# AssignEnvironmentCondition().AssignCondition(solidification_model_part)  
    
    BiphasicFillingUtilities().ComputeNodalVolume(solidification_model_part) 
    
    solidification_model_part.ProcessInfo = fluid_model_part.ProcessInfo
    solidification_them_solver = convdiff_phasechange_solver.Solver(
        solidification_model_part, domain_size, my_settings)
    #solidification_them_solver.MaxNewtonRapshonIterations = 15
    #solidification_them_solver.MaxLineSearchIterations = 10
    solidification_them_solver.use_linesearch_in_step1=True
    solidification_them_solver.skip_stage0 = True
    solidification_them_solver.stage1_max_iterations = 5
    solidification_them_solver.Initialize()

    for node in fluid_model_part.Nodes:
        node.Free(TEMPERATURE)
        node.SetSolutionStepValue(VELOCITY_X, 0.0)
        node.SetSolutionStepValue(VELOCITY_Y, 0.0)
        node.SetSolutionStepValue(VELOCITY_Z, 0.0)
        node.SetSolutionStepValue(DISTANCE, -1.0)
        node.SetSolutionStepValue(IS_VISITED, 0.0)
    import os
    if(os.path.isfile(input_file_name + ".Tmp")):
        filling_temp = open(input_file_name + ".Tmp", 'r')
        for line in filling_temp:
            ll = line.strip()
            clm = ll.split()
            nd_id = clm[0]
            nd_temp = clm[1]
            fluid_model_part.Nodes[int(nd_id)].SetSolutionStepValue(
                TEMPERATURE, float(nd_temp))
        filling_temp.close()

    cooled_precent = 0.0
    solidification_model_part.ProcessInfo.SetValue(IS_SOLIDIFIED, 0)
    solidification_percent = 0.02
    max_cooling_delta_temp = 30.0
    solidification_time = EstimateTimeStep3D().EstimateSolidificationTime(
        solidification_model_part)
    print (solidification_time)

    dt_min_ss = solidification_time / 500.0
    dt_max_ss = solidification_time / 100.0
    Dt = 0.1 * dt_min_ss
    stop_temp = ProjectParameters.AMBIENT_TEMPERATURE

    shrinkage_porosity_calculation = ShrinkagePorosityCalculationProcess(
        fluid_model_part)

    while(cooled_precent <= 99.5 and time < 2.0 * solidification_time):

        time = time + Dt
        solidification_model_part.CloneTimeStep(time)
        elapsed_time = time_measure.time() - start_time
        print ("Solidifaction_DT= ", Dt, "max_dt= ",
               dt_max_ss, "Total time= ", elapsed_time)
        solidification_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.0)

        solidification_them_solver.Solve()
        
        
        #for node in fluid_model_part.Nodes:
            ##node.SetSolutionStepValue(DENSITY, 0, initial_density)
            #print(node.GetSolutionStepValue(TEMPERATURE))
            #print(node.GetSolutionStepValue(DENSITY))
        
        shrinkage_porosity_calculation.Execute()
        Dt = EstimateTimeStep3D().ComputeSolidificationCoolingDt(
            solidification_model_part,  solidification_percent,  max_cooling_delta_temp,  dt_min_ss,  dt_max_ss)
        cooled_precent = EstimateTimeStep3D().CheckStopTemperature(
            solidification_model_part, stop_temp)

        if(cooled_precent >= next_screen_output):
            print("*********************************************")
            print("Solidification-Cooling %.0f" % cooled_precent, "%\t")
            print("*********************************************")
            sys.stdout.flush()
            next_screen_output += screen_output_dt

        if(time >= next_output_time):
            verifyp.Execute()
            copy_solidification_variables_process.Execute()
            gid_io.WriteNodalResults(
                TEMPERATURES, fluid_model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(
                SOLIDIF_TIME, fluid_model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(
                SOLIDIF_MODULUS, fluid_model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(
                DENSITY, fluid_model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(
                HTC, fluid_model_part.Nodes, time, 0)
           
            gid_io.WriteNodalResults(
                SOLID_FRACTION, fluid_model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(
                MACRO_POROSITY, fluid_model_part.Nodes, time, 0)
            if(write_US_units):
                gid_io.WriteNodalResults(
                    SOLIDIF_MODULUS_US, fluid_model_part.Nodes, time, 0)
                gid_io.WriteNodalResults(
                    TEMPERATURES_US, fluid_model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(
                TEMPERATURE, fluid_model_part.Nodes, time, 0)
            # gid_io.WriteNodalResults(SOLIDFRACTION,fluid_model_part.Nodes,time,0)
            last_print = time
            gid_io.Flush()
            sys.stdout.flush()

            out = 0
            next_output_time += output_dt
        out = out + 1
        step = step + 1
    # Special treatment for the SHRINKAGE_POROSITY
    non_zero_shrinkage_nodes = NodesArray()
    BiphasicFillingUtilities().MacroPorosityToShrinkageComputation(
        fluid_model_part, non_zero_shrinkage_nodes, 10)
    copy_solidification_variables_process.Execute()
    gid_io.WriteNodalResults(
        SHRINKAGE_POROSITY, non_zero_shrinkage_nodes, time, 0)
    if(write_US_units):
        gid_io.WriteNodalResults(
            SHRINKAGE_POROSITY_US, non_zero_shrinkage_nodes, time, 0)

    # write the last step results
    if(last_print != time):
# copy_solidification_variables_process.Execute()
        gid_io.WriteNodalResults(TEMPERATURES, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(SOLIDIF_TIME, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(
            SOLIDIF_MODULUS, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(
            SOLID_FRACTION, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(
            MACRO_POROSITY, fluid_model_part.Nodes, time, 0)
        if(write_US_units):
            gid_io.WriteNodalResults(
                SOLIDIF_MODULUS_US, fluid_model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(
                TEMPERATURES_US, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(TEMPERATURE, fluid_model_part.Nodes, time, 0)
        # gid_io.WriteNodalResults(SOLIDFRACTION,fluid_model_part.Nodes,time,0)

        gid_io.FinalizeResults()
        gid_io.Flush()
        sys.stdout.flush()
    gid_io.FinalizeResults()
    print ("Click2cast running state: End SOLIDIFICATION-COOLING solver")
    sys.stdout.flush()


print("Click2cast running state: End KRATOS")
sys.stdout.flush()
