from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

kratos_path = '../../../..'
kratos_benchmarking_path = '../../../../benchmarking'  # kratos_root/benchmarking
import sys
sys.path.append(kratos_path)
sys.path.append(kratos_benchmarking_path)

from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FSIApplication import *

import benchmarking
import NonConformant_OneSideMap

def BenchmarkCheck(o_model_part, d_model_part):
    # Find some nodes in both meshes and compare their vaules
    for node in o_model_part.Nodes:
        if (node.X < 0.51) and (node.X > 0.49) and (node.Y < 0.51) and (node.Y > 0.49) and (node.Z < 0.51) and (node.Z > 0.49):
            press_1o = node.GetSolutionStepValue(PRESSURE, 0)
            velX_1o = node.GetSolutionStepValue(VELOCITY_X, 0)
            velY_1o = node.GetSolutionStepValue(VELOCITY_Y, 0)
            velZ_1o = node.GetSolutionStepValue(VELOCITY_Z, 0)

    for node in d_model_part.Nodes:
        if (node.X < 0.51) and (node.X > 0.49) and (node.Y < 0.51) and (node.Y > 0.49) and (node.Z < 0.51) and (node.Z > 0.49):
            press_1d = node.GetSolutionStepValue(PRESSURE, 0)
            velX_1d = node.GetSolutionStepValue(VELOCITY_X, 0)
            velY_1d = node.GetSolutionStepValue(VELOCITY_Y, 0)
            velZ_1d = node.GetSolutionStepValue(VELOCITY_Z, 0)

    press_1 = press_1d - press_1o
    
    velX_1 = velX_1o - velX_1d
    velY_1 = velY_1o - velY_1d
    velZ_1 = velZ_1o - velZ_1d

    benchmarking.Output(press_1, "Difference in pressure in test point 1")
    benchmarking.Output(velX_1, "Difference in x velocity in test point 1")
    benchmarking.Output(velY_1, "Difference in y velocity in test point 1")
    benchmarking.Output(velZ_1, "Difference in z velocity in test point 1")

# Json format solvers settings
fluid_parameter_file = open("ProjectParametersFluid.json",'r')
ProjectParametersFluid = Parameters( fluid_parameter_file.read())
solid_parameter_file = open("ProjectParametersSolid.json",'r')
ProjectParametersSolid = Parameters( solid_parameter_file.read())

# Defining a model part for the fluid and one for the structure
origin_model_part = ModelPart("structure_part")
destination_model_part = ModelPart("fluid_part")

# Set the domain size (2D test)
origin_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, 3)
destination_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, 3)

# Create the fluid and solid solvers
fluid_solver_module = __import__(ProjectParametersFluid["solver_settings"]["solver_type"].GetString())
fluid_solver = fluid_solver_module.CreateSolver(destination_model_part,ProjectParametersFluid["solver_settings"])
solid_solver_module = __import__(ProjectParametersSolid["solver_settings"]["solver_type"].GetString())
solid_solver = solid_solver_module.CreateSolver(origin_model_part,ProjectParametersSolid["solver_settings"])

# Variables addition
fluid_solver.AddVariables()
fluid_solver.main_model_part.AddNodalSolutionStepVariable(FORCE) # Specific addition of FORCE variable to check the nodal fluxes
solid_solver.AddVariables()
NonConformant_OneSideMap.AddVariables(destination_model_part, origin_model_part)

# Import model_parts
fluid_solver.ImportModelPart()
solid_solver.ImportModelPart()

# Add DOFs
fluid_solver.AddDofs()
solid_solver.AddDofs()

# Initialize solvers
fluid_solver.Initialize()
solid_solver.Initialize()

# Initialize GidIO
input_file = "meshtest_3D"
gid_mode = GiDPostMode.GiD_PostAscii
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(input_file, gid_mode, multifile, deformed_mesh_flag, write_conditions)

for node in origin_model_part.GetSubModelPart("Structure_interface").Nodes:
    node.SetSolutionStepValue(VELOCITY_X, 0, 1.0)
    node.SetSolutionStepValue(VELOCITY_Y, 0, 2.0)
    node.SetSolutionStepValue(VELOCITY_Z, 0, 3.0)
    node.Set(INTERFACE)
print("***** Origin mesh data assigned to nodes *****")

for node in destination_model_part.GetSubModelPart("Fluid_interface").Nodes:
    node.SetSolutionStepValue(PRESSURE, 0, 4.0)
    node.SetSolutionStepValue(FORCE_X, 0, 5.0)
    node.SetSolutionStepValue(FORCE_Y, 0, 6.0)
    node.SetSolutionStepValue(FORCE_Z, 0, 7.0)
    node.Set(INTERFACE)
print("***** Destination mesh data assigned to nodes *****")

# Print original mesh
gid_io.InitializeMesh(0)
gid_io.WriteMesh(origin_model_part.GetMesh())
gid_io.FinalizeMesh()

gid_io.InitializeResults(0, origin_model_part.GetMesh())
gid_io.WriteNodalResults(PRESSURE, origin_model_part.Nodes, 0, 0)
gid_io.WriteNodalResults(VELOCITY, origin_model_part.Nodes, 0, 0)
gid_io.FinalizeResults()

print(destination_model_part)
print(origin_model_part)

search_radius_factor = 1.0
mapper_max_iteration = 200
mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(destination_model_part, 
                                                           origin_model_part, 
                                                           search_radius_factor,
                                                           mapper_max_iteration)
print("***** NonConformant mapper created *****")

print("***** Attempting transfer *****")
# Transfer pressure to destination mesh
mapper.FluidToStructure_ScalarMap(PRESSURE, PRESSURE, True)
mapper.FluidToStructure_VectorMap(FORCE, POINT_LOAD, True, True)
mapper.StructureToFluid_VectorMap(VELOCITY, VELOCITY, True, False)
print("***** Information transferred *****")

if (benchmarking.InBuildReferenceMode()):
    # Output expected values
    BenchmarkCheck(origin_model_part, destination_model_part)
else:
    # Output calculated values
    BenchmarkCheck(origin_model_part, destination_model_part)

# Print destination mesh
gid_io.InitializeMesh(1)
gid_io.WriteMesh(destination_model_part.GetMesh())
gid_io.FinalizeMesh()

gid_io.InitializeResults(1, destination_model_part.GetMesh())
gid_io.WriteNodalResults(PRESSURE, destination_model_part.Nodes, 0, 0)
gid_io.WriteNodalResults(FORCE, destination_model_part.Nodes, 0, 0)
gid_io.WriteNodalResults(VAUX_EQ_TRACTION, destination_model_part.Nodes, 0, 0)
gid_io.WriteNodalResults(NODAL_MAUX, destination_model_part.Nodes, 0, 0)
gid_io.WriteNodalResults(VELOCITY, destination_model_part.Nodes, 0, 0)
gid_io.FinalizeResults()

# Print (check) origin mesh (note that PRESSURE is mapped from the destiny mesh to the origin one)
gid_io.InitializeMesh(2)
gid_io.WriteMesh(origin_model_part.GetMesh())
gid_io.FinalizeMesh()

gid_io.InitializeResults(2, origin_model_part.GetMesh())
gid_io.WriteNodalResults(PRESSURE, origin_model_part.Nodes, 0, 0)
gid_io.WriteNodalResults(POINT_LOAD, origin_model_part.Nodes, 0, 0)
gid_io.WriteNodalResults(VAUX_EQ_TRACTION, origin_model_part.Nodes, 0, 0)
gid_io.WriteNodalResults(NODAL_MAUX, origin_model_part.Nodes, 0, 0)
gid_io.WriteNodalResults(VELOCITY, origin_model_part.Nodes, 0, 0)
gid_io.FinalizeResults()

# Interface fluxes equilibrium assessment
sum_fx_solid = 0.0
sum_fy_solid = 0.0
sum_fz_solid = 0.0

for node in origin_model_part.GetSubModelPart("Structure_interface").Nodes:
    f_solid = solid_solver.main_model_part.Nodes[node.Id].GetSolutionStepValue(POINT_LOAD)
    sum_fx_solid += f_solid[0]
    sum_fy_solid += f_solid[1]
    sum_fz_solid += f_solid[2]
    
sum_fx_fluid = 0.0
sum_fy_fluid = 0.0
sum_fz_fluid = 0.0

for node in destination_model_part.GetSubModelPart("Fluid_interface").Nodes:
    f_fluid = fluid_solver.main_model_part.Nodes[node.Id].GetSolutionStepValue(FORCE)
    sum_fx_fluid += f_fluid[0]
    sum_fy_fluid += f_fluid[1]
    sum_fz_fluid += f_fluid[2]
    
print("")
print("### INTERFACE FLUXES EQUILIBRIUM RESULTS ###")
print("Sum. Fx fluid: ", sum_fx_fluid, " Sum. Fx solid: ", sum_fx_solid)
print("Sum. Fy fluid: ", sum_fy_fluid, " Sum. Fy solid: ", sum_fy_solid)
print("Sum. Fz fluid: ", sum_fz_fluid, " Sum. Fz solid: ", sum_fz_solid)
print("Relative err. Fx: ",abs(100*(sum_fx_fluid-sum_fx_solid)/(sum_fx_fluid+10e-8)))
print("Relative err. Fy: ",abs(100*(sum_fy_fluid-sum_fy_solid)/(sum_fy_fluid+10e-8)))
print("Relative err. Fz: ",abs(100*(sum_fz_fluid-sum_fz_solid)/(sum_fz_fluid+10e-8)))
