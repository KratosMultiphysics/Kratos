import mesh_only_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = mesh_only_var.domain_size

kratos_path = '../../../..'
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
import sys
sys.path.append(kratos_path)
sys.path.append(kratos_benchmarking_path)

from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FSIApplication import *

##################################################################
##################################################################

import benchmarking

def BenchmarkCheck(o_model_part, d_model_part):
    # Find some nodes in both meshes and compare their vaules
    for node in o_model_part.Nodes:
        if (node.X == -5) and (node.Y == -5) and (node.Z == -5):
            press_1o = node.GetSolutionStepValue(PRESSURE,0)
            velX_1o = node.GetSolutionStepValue(VELOCITY_X,0)
            velY_1o = node.GetSolutionStepValue(VELOCITY_Y,0)
            velZ_1o = node.GetSolutionStepValue(VELOCITY_Z,0)
        elif (node.X == -5) and (node.Y == 5) and (node.Z == -5):
            press_2o = node.GetSolutionStepValue(PRESSURE,0)
            velX_2o = node.GetSolutionStepValue(VELOCITY_X,0)
            velY_2o = node.GetSolutionStepValue(VELOCITY_Y,0)
            velZ_2o = node.GetSolutionStepValue(VELOCITY_Z,0)
        elif (node.X == 5) and (node.Y == 5) and (node.Z == 5):
            press_3o = node.GetSolutionStepValue(PRESSURE,0)
            velX_3o = node.GetSolutionStepValue(VELOCITY_X,0)
            velY_3o = node.GetSolutionStepValue(VELOCITY_Y,0)
            velZ_3o = node.GetSolutionStepValue(VELOCITY_Z,0)
        elif (node.X == 5) and (node.Y == -5) and (node.Z == 5):
            press_4o = node.GetSolutionStepValue(PRESSURE,0)
            velX_4o = node.GetSolutionStepValue(VELOCITY_X,0)
            velY_4o = node.GetSolutionStepValue(VELOCITY_Y,0)
            velZ_4o = node.GetSolutionStepValue(VELOCITY_Z,0)

    for node in d_model_part.Nodes:
        if (node.X == -5) and (node.Y == -5) and (node.Z == -5):
            press_1d = node.GetSolutionStepValue(PRESSURE,0)
            velX_1d = node.GetSolutionStepValue(VELOCITY_X,0)
            velY_1d = node.GetSolutionStepValue(VELOCITY_Y,0)
            velZ_1d = node.GetSolutionStepValue(VELOCITY_Z,0)
        elif (node.X == -5) and (node.Y == 5) and (node.Z == -5):
            press_2d = node.GetSolutionStepValue(PRESSURE,0)
            velX_2d = node.GetSolutionStepValue(VELOCITY_X,0)
            velY_2d = node.GetSolutionStepValue(VELOCITY_Y,0)
            velZ_2d = node.GetSolutionStepValue(VELOCITY_Z,0)
        elif (node.X == 5) and (node.Y == 5) and (node.Z == 5):
            press_3d = node.GetSolutionStepValue(PRESSURE,0)
            velX_3d = node.GetSolutionStepValue(VELOCITY_X,0)
            velY_3d = node.GetSolutionStepValue(VELOCITY_Y,0)
            velZ_3d = node.GetSolutionStepValue(VELOCITY_Z,0)
        elif (node.X == 5) and (node.Y == -5) and (node.Z == 5):
            press_4d = node.GetSolutionStepValue(PRESSURE,0)
            velX_4d = node.GetSolutionStepValue(VELOCITY_X,0)
            velY_4d = node.GetSolutionStepValue(VELOCITY_Y,0)
            velZ_4d = node.GetSolutionStepValue(VELOCITY_Z,0)

    press_1 = press_1d - press_1o
    press_2 = press_2d - press_2o
    press_3 = press_3d - press_3o
    press_4 = press_4d - press_4o

    velX_1 = velX_1o - velX_1d
    velX_2 = velX_2o - velX_2d
    velX_3 = velX_3o - velX_3d
    velX_4 = velX_4o - velX_4d

    velY_1 = velY_1o - velY_1d
    velY_2 = velY_2o - velY_2d
    velY_3 = velY_3o - velY_3d
    velY_4 = velY_4o - velY_4d

    velZ_1 = velZ_1o - velZ_1d
    velZ_2 = velZ_2o - velZ_2d
    velZ_3 = velZ_3o - velZ_3d
    velZ_4 = velZ_4o - velZ_4d

    benchmarking.Output(press_1, "Difference in pressure in test point 1")
    benchmarking.Output(velX_1, "Difference in x velocity in test point 1")
    benchmarking.Output(velY_1, "Difference in y velocity in test point 1")
    benchmarking.Output(velZ_1, "Difference in z velocity in test point 1")

    benchmarking.Output(press_2, "Difference in pressure in test point 2")
    benchmarking.Output(velX_2, "Difference in x velocity in test point 2")
    benchmarking.Output(velY_2, "Difference in y velocity in test point 2")
    benchmarking.Output(velZ_2, "Difference in z velocity in test point 2")

    benchmarking.Output(press_3, "Difference in pressure in test point 3")
    benchmarking.Output(velX_3, "Difference in x velocity in test point 3")
    benchmarking.Output(velY_3, "Difference in y velocity in test point 3")
    benchmarking.Output(velZ_3, "Difference in z velocity in test point 3")
    
    benchmarking.Output(press_4, "Difference in pressure in test point 4")
    benchmarking.Output(velX_4, "Difference in x velocity in test point 4")
    benchmarking.Output(velY_4, "Difference in y velocity in test point 4")
    benchmarking.Output(velZ_4, "Difference in z velocity in test point 4")
    
#defining a model part for the fluid and one for the structure
origin_model_part = ModelPart("structure_part")
destination_model_part = ModelPart("fluid_part")

#############################################
##importing the solvers needed
import fractional_step_solver
import NonConformant_OneSideMap
import structural_solver_dynamic

fractional_step_solver.AddVariables(destination_model_part)
NonConformant_OneSideMap.AddVariables(destination_model_part,origin_model_part)
structural_solver_dynamic.AddVariables(origin_model_part)

#introducing input file name
input_file = mesh_only_var.problem_name
destination_mesh = "destination2"

#reading the original part
gid_mode = GiDPostMode.GiD_PostBinary
##gid_mode = GiDPostMode.GiD_PostAscii
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(input_file,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_origin = ModelPartIO(input_file)
model_part_io_origin.ReadModelPart(origin_model_part)
print "***** Origin mesh read *****"

#reading destination mesh
model_part_io_destination = ModelPartIO(destination_mesh)
model_part_io_destination.ReadModelPart(destination_model_part)
print "***** Destination mesh read *****"

# Assign a variable pressure to the mesh and mark nodes as interface
for node in origin_model_part.Nodes:
    node.SetSolutionStepValue(PRESSURE,0,node.Y+node.Z)
    node.SetSolutionStepValue(IS_INTERFACE,0,1.0)
print "***** Pressure assigned to nodes *****"

for node in destination_model_part.Nodes:
    node.SetSolutionStepValue(VELOCITY_X,0,-node.X-node.Z)
    node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
    node.SetSolutionStepValue(VELOCITY_Z,0,node.X+node.Y)
    node.SetSolutionStepValue(IS_INTERFACE,0,1.0)
print "***** Velocity assigned to nodes *****"

# Print original mesh
gid_io.InitializeMesh( 0 )
gid_io.WriteMesh( origin_model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults( 0 ,origin_model_part.GetMesh())
gid_io.WriteNodalResults(PRESSURE,origin_model_part.Nodes,0,0)
gid_io.WriteNodalResults(VELOCITY,origin_model_part.Nodes,0,0)
gid_io.WriteNodalResults(IS_INTERFACE,origin_model_part.Nodes,0,0)
gid_io.FinalizeResults()

print "Attempting transfer"
# Transfer pressure to destination mesh
mapper=NonConformant_OneSideMap.\
        NonConformant_OneSideMap(destination_model_part,origin_model_part,\
                                 1.0,15)
print "***** mapper created *****"

mapper.StructureToFluid_ScalarMap(PRESSURE,PRESSURE)
mapper.FluidToStructure_VectorMap(VELOCITY,VELOCITY)
print "***** information transferred *****"

if (benchmarking.InBuildReferenceMode()):
    # Output expected values
    BenchmarkCheck(origin_model_part,destination_model_part)
else:
    # Output calculated values
    BenchmarkCheck(origin_model_part,destination_model_part)

# Print destination mesh
gid_io.InitializeMesh( 1 )
gid_io.WriteMesh( destination_model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults( 1 ,destination_model_part.GetMesh())
gid_io.WriteNodalResults(PRESSURE,destination_model_part.Nodes,0,0)
gid_io.WriteNodalResults(VELOCITY,destination_model_part.Nodes,0,0)
gid_io.WriteNodalResults(IS_INTERFACE,destination_model_part.Nodes,0,0)
gid_io.FinalizeResults()

        
# Print (check) origin mesh
gid_io.InitializeMesh( 2 )
gid_io.WriteMesh( origin_model_part.GetMesh() )
gid_io.FinalizeMesh()

gid_io.InitializeResults( 2 ,origin_model_part.GetMesh())
gid_io.WriteNodalResults(PRESSURE,origin_model_part.Nodes,0,0)
gid_io.WriteNodalResults(VELOCITY,origin_model_part.Nodes,0,0)
gid_io.WriteNodalResults(IS_INTERFACE,origin_model_part.Nodes,0,0)
gid_io.FinalizeResults()
