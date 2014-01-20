from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
domain_size = 3

from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *

import math
import time

# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")

proc_info = fluid_model_part.ProcessInfo

# defining a model part for the fluid and one for the structure
fluid_model_part.AddNodalSolutionStepVariable(RADIUS)
fluid_model_part.AddNodalSolutionStepVariable(NODAL_AREA)
fluid_model_part.AddNodalSolutionStepVariable(NODAL_MASS)
fluid_model_part.AddNodalSolutionStepVariable(RHS)
fluid_model_part.AddNodalSolutionStepVariable(WORK)
fluid_model_part.AddNodalSolutionStepVariable(FLOW)
fluid_model_part.AddNodalSolutionStepVariable(VELOCITY)
fluid_model_part.AddNodalSolutionStepVariable(THICKNESS)
fluid_model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS)
fluid_model_part.AddNodalSolutionStepVariable(POISSON_RATIO)
fluid_model_part.AddNodalSolutionStepVariable(DENSITY)
fluid_model_part.AddNodalSolutionStepVariable(TERMINAL_RESISTANCE)
fluid_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)


# introducing input file name
# input_file_name = "HeartModel16arteries3parts_Right_Balanced_Dominant_coupling"
input_file_name = "alastruey_new_coronaries"
# input_file_name = "test_alas"

# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(2)

# import removal_tool
# removal_tool.DoRemoval(model_part1D)

# for element in fluid_model_part.Elements:
 # properties = element.Properties
 # print properties.Id
 # E = properties[YOUNG_MODULUS]
 # h = properties[THICKNESS]
 # nu= properties[POISSON_RATIO]
 # den=properties[DENSITY]
 # r=properties[RADIUS]
 # area_inicial=r*r*3.141617
 # for node in element.GetNodes():
   # node.SetSolutionStepValue(YOUNG_MODULUS,E)
   # node.SetSolutionStepValue(THICKNESS,h)
   # node.SetSolutionStepValue(POISSON_RATIO,nu)
   # node.SetSolutionStepValue(DENSITY,den)
   # node.SetSolutionStepValue(RADIUS,r)
   # node.SetSolutionStepValue(NODAL_AREA,area_inicial)

# for node in fluid_model_part.Nodes:
    # print node.Id
    # AREA= node.GetSolutionStepValue(NODAL_AREA)
    # h= node.GetSolutionStepValue(THICKNESS)
    # nu= node.GetSolutionStepValue(POISSON_RATIO)
    # raidus= node.GetSolutionStepValue(RADIUS)
    # print AREA
    # print E
    # print h
    # print nu
    # raw_input()


# settings to be changed
Dt = 1e-5
full_Dt = Dt
initial_Dt = Dt  # 0.001 * full_Dt #0.05 #0.01
cardiac_cycle = 3.0
time_cardiac_cycle = 0.83
final_time = cardiac_cycle * time_cardiac_cycle
output_step = 100

out = 1

# mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh(fluid_model_part.GetMesh())
gid_io.FinalizeMesh()
gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())

integrator = ArteryTimeIntegrator()
integrator.Initialize(fluid_model_part)
inletconditiontable = fluid_model_part.GetTable(3)

# Aqui estoy inizializando mis variables
q_inicial = inletconditiontable.GetValue(0)

for node in fluid_model_part.Nodes:
    if(node.IsFixed(FLOW) == False):
    # node.SetSolutionStepValue(FLOW,0,0.0008)
        node.SetSolutionStepValue(FLOW, 0, q_inicial)
        # node.SetSolutionStepValue(NODAL_AREA,0,0.0006157521594)

    else:
        # node.SetSolutionStepValue(FLOW,0,0.0008)
        node.SetSolutionStepValue(FLOW, 0, q_inicial)
        # properties = element.Properties
        # r=properties[RADIUS]
        # area_inicial=r*r*3.141617
        # node.SetSolutionStepValue(NODAL_AREA,0,0.0006157521594)
        # node.SetSolutionStepValue(NODAL_AREA,area_inicial)


# for node in fluid_model_part.Nodes:
    # area = node.GetSolutionStepValue(NODAL_AREA)
    # flow=node.GetSolutionStepValue(FLOW)
    # radius=node.GetSolutionStepValue(NODAL_AREA)
    # print "area=", area
    # print "FLOW", flow
time = 0.0
total_time = 0.0
step = 0
while(total_time < final_time):
        # print "--------------------------01062013---------------------------------------------------------------"
        # raw_input()

#    if(step < 5):
#        Dt = initial_Dt
#    else:
#        Dt = full_Dt
        # Dt = integrator.EstimateDeltaTime(fluid_model_part, 0.4)
        # Dt = 1e-4
        # print "Delta_time::::>>", Dt

        # time = time + Dt
    print("time=", time)
    # Q=inletconditiontable.GetValue(time)
    # Q=0.0
    # print "Caudal que estoy imponiendo=", Q
    # node.SetSolutionStepValue(FLOW,Q)
    # if(time > final_time / 3.0):
    for node in fluid_model_part.Nodes:
        if(node.IsFixed(FLOW)):
            # puls= math.sin(2.00*math.pi*time/(1.))
            # if(time > 0.0005 and time < 0.5):
              # puls = 1.00
            # else:
              # puls = 0.0
            # velocity = 0.0008*puls
            # Q=velocity*area
            # Q=0.0008*puls
            # Q=0.0
            Q = inletconditiontable.GetValue(time)
            node.SetSolutionStepValue(FLOW, 0, Q)
            # print "Caudal que estoy imponiendo=", Q, "tiempo-->", time, "nodo-->", node.Id
            # cau = node.GetSolutionStepValue(FLOW,1)
            # print "caudal===>", cau, "nodo", node
            if(step > 1):
                Q = 2 * Q - node.GetSolutionStepValue(FLOW, 1)
                # print "Caudal que estoy imponiendo(RIEMMAANN)=", Q, "tiempo-->", time, "nodo-->", node
                node.SetSolutionStepValue(FLOW, 0, Q)

    # print time
    # for node in fluid_model_part.Nodes:
        # area = node.GetSolutionStepValue(NODAL_AREA)
        # print "area===>", area, "nodo", node
    # raw_input()
    # for node in fluid_model_part.Nodes:
        # area = node.GetSolutionStepValue(NODAL_AREA)
        # flow=node.GetSolutionStepValue(FLOW)
        # print "area=", area, "FLOW", flow, "nodo", node
    fluid_model_part.CloneTimeStep(total_time)
    integrator.SolveStep(fluid_model_part)

    # for node in fluid_model_part.Nodes:
        # area = node.GetSolutionStepValue(NODAL_AREA)
        # print "area===>", area, "nodo", node.Id
    # raw_input()
    # raw_input()

    # for node in fluid_model_part.Nodes:
        # area = node.GetSolutionStepValue(NODAL_AREA)
        # flow=node.GetSolutionStepValue(FLOW)
        # radius=node.GetSolutionStepValue(NODAL_AREA)
        # print "area=", area
        # print "FLOW", flow

    # raw_input()

    if(out == output_step):
    # para escribir por Dt
    # if(output_time <= out):
        gid_io.WriteNodalResults(NODAL_AREA, fluid_model_part.Nodes, total_time, 0)
        gid_io.WriteNodalResults(NODAL_MASS, fluid_model_part.Nodes, total_time, 0)
        gid_io.WriteNodalResults(FLOW, fluid_model_part.Nodes, total_time, 0)
        gid_io.WriteNodalResults(RADIUS, fluid_model_part.Nodes, total_time, 0)
        gid_io.WriteNodalResults(RHS, fluid_model_part.Nodes, total_time, 0)
        gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, total_time, 0)
        gid_io.WriteNodalResults(WORK, fluid_model_part.Nodes, total_time, 0)
        # gid_io.WriteNodalResults(NODAL_MASS,fluid_model_part.Nodes,total_time,0)

        out = 0

    # para escribir en funcion de dt
    # out = out + Dt
    # print "final_time=",time_cardiac_cycle

    if(time > time_cardiac_cycle):
        time = time - time_cardiac_cycle

    Dt = 1e-4
    # print "Delta_time::::>>", Dt
    time = time + Dt
    total_time = total_time + Dt
    out = out + 1
    step = step + 1

gid_io.FinalizeResults()
