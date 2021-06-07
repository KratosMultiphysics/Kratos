# time control starts
import time as timer
print(timer.ctime())

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.TopologyOptimizationApplication import *
from KratosMultiphysics.LinearSolversApplication import *

# For GID output
from KratosMultiphysics.gid_output import GiDOutput

print("checkpoint 0\n"+"-"*50)

# GiD output settings
nodal_results=[]
gauss_points_results=["X_PHYS"]
VolumeOutput = True
GiDPostMode = "Binary"
GiDWriteMeshFlag = False
GiDWriteConditionsFlag = True
GiDWriteParticlesFlag = False
GiDMultiFileFlag = "Single"
       

print("checkpoint 1\n"+"-"*50)
# Read optimized model part from restart file
model = Model()
optimized_model_part = model.CreateModelPart("optimized_model_part")
optimized_model_part.AddNodalSolutionStepVariable(NORMAL)
restart_file_name = "Small_Cantilever_Restart_File_6"
model_part_io = ModelPartIO(restart_file_name)
model_part_io.ReadModelPart(optimized_model_part)

print("checkpoint 2\n"+"-"*50)
# Extract volume
threshold = 0.2
model = Model()
extracted_volume_model_part = model.CreateModelPart("extracted_volume_model_part")
TopologyExtractorUtilities().ExtractVolumeMesh(optimized_model_part, 
                                               threshold, 
                                               extracted_volume_model_part)


print("checkpoint 3\n"+"-"*50)
# Write extracted volume
gid_io = GiDOutput("extracted_volume_model_part", VolumeOutput, GiDPostMode, 
                   GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
gid_io.initialize_results(extracted_volume_model_part)
# gid_io.write_results(1, extracted_volume_model_part, nodal_results, gauss_points_results)
# gid_io.finalize_results()

# print("checkpoint 4\n"+"-"*50)
# # Extract surface
# model = Model()
# extracted_surface_model_part = model.CreateModelPart("extracted_surface_model_part")
# TopologyExtractorUtilities().ExtractSurfaceMesh(extracted_volume_model_part, extracted_surface_model_part)

# print("checkpoint 5\n"+"-"*50)
# # Write extracted surface
# gid_io_2 = GiDOutput("extracted_surface_model_part", VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
# gid_io_2.initialize_results(extracted_surface_model_part)
# gid_io_2.write_results(1, extracted_surface_model_part, nodal_results, gauss_points_results)
# gid_io_2.finalize_results()

# print("checkpoint 6\n"+"-"*50)
# # Smooth extracted surface
# smoothing_relaxation_factor = 1
# smoothing_iterations = 5
# model = Model()
# smoothed_surface_model_part = model.CreateModelPart("extracted_surface_model_part")
# TopologySmoothingUtilities().SmoothMesh(extracted_surface_model_part, smoothing_relaxation_factor, smoothing_iterations)

# # Write smoothed mesh
# gid_io_3 = GiDOutput("smoothed_surface_model_part", VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
# gid_io_3.initialize_results(extracted_surface_model_part)
# gid_io_3.write_results(1, extracted_surface_model_part, nodal_results, gauss_points_results)
# gid_io_3.finalize_results()

# # Write stl of extracted surface
# IOUtilities().WriteSurfaceAsSTLFile("smoothed_design_21.stl",extracted_surface_model_part)
