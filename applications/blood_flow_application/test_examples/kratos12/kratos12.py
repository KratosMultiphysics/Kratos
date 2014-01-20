from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
domain_size = 3

from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *

import math

# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")

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
input_file_name = "kratos12"

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

for element in fluid_model_part.Elements:
    properties = element.Properties
    print(properties.Id)
    E = properties[YOUNG_MODULUS]
    h = properties[THICKNESS]
    nu = properties[POISSON_RATIO]
    den = properties[DENSITY]
    r = properties[RADIUS]
    area_inicial = r * r * 3.141617
    for node in element.GetNodes():
        node.SetSolutionStepValue(YOUNG_MODULUS, E)
        node.SetSolutionStepValue(THICKNESS, h)
        node.SetSolutionStepValue(POISSON_RATIO, nu)
        node.SetSolutionStepValue(DENSITY, den)
        node.SetSolutionStepValue(RADIUS, r)
        node.SetSolutionStepValue(NODAL_AREA, area_inicial)


for node in fluid_model_part.Nodes:
    if(node.IsFixed(FLOW) == False):
        node.SetSolutionStepValue(FLOW, 0, 0.0008)
        node.SetSolutionStepValue(NODAL_AREA, 0, 0.0006157521594)

    else:
        node.SetSolutionStepValue(FLOW, 0, 0.0008)
        node.SetSolutionStepValue(NODAL_AREA, 0, 0.0006157521594)

# for node in fluid_model_part.Nodes:
#    node.SetSolutionStepValue(FLOW,0,0.0)

# settings to be changed
Dt = 4e-5
full_Dt = Dt
initial_Dt = Dt  # 0.001 * full_Dt #0.05 #0.01
final_time = 5
output_step = 5

out = 1

# mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh(fluid_model_part.GetMesh())
gid_io.FinalizeMesh()
gid_io.InitializeResults(mesh_name, (fluid_model_part).GetMesh())

integrator = ArteryTimeIntegrator()
integrator.Initialize(fluid_model_part)


time = 0.0
step = 0
while(time < final_time):
    print("--------------------------------------------------------------------------------------------------------")

    if(step < 5):
        Dt = initial_Dt
    else:
        Dt = full_Dt

    # Dt = integrator.EstimateDeltaTime(fluid_model_part, 0.4)
    Dt = 4e-5
    print("Dt", Dt)

    time = time + Dt
    print("time=", time)
#    if(time > final_time / 3.0):
    for node in fluid_model_part.Nodes:
        puls = math.sin(2.00 * math.pi * time / (1.))
#      if(time > 0.25 and time < 0.5):
#        puls = 1.00
#      else:
        if puls < 0.00:
            puls = 0.0
        if(node.IsFixed(FLOW)):
            # velocity = 0.0008*puls
            # area = node.GetSolutionStepValue(NODAL_AREA)
            # Q=velocity*area
            Q = 0.0008 * puls
            node.SetSolutionStepValue(FLOW, 0, Q)

    print(time)
    fluid_model_part.CloneTimeStep(time)

    integrator.SolveStep(fluid_model_part)

    if(out == output_step):
        gid_io.WriteNodalResults(NODAL_AREA, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(NODAL_MASS, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(FLOW, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(RADIUS, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(RHS, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(WORK, fluid_model_part.Nodes, time, 0)
        # gid_io.WriteNodalResults(NODAL_MASS,fluid_model_part.Nodes,time,0)

        out = 0

    out = out + 1
    step = step + 1

gid_io.FinalizeResults()
