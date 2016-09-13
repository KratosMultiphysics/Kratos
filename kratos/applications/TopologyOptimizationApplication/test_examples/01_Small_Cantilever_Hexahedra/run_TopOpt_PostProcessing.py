from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#### TIME MONITORING START ####

# time control starts
import time as timer
print(timer.ctime())

#### TIME MONITORING END ####

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.TopologyOptimizationApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

# For GID output
from gid_output import GiDOutput

######################################################################################
######################################################################################
######################################################################################

# GiD output settings
nodal_results=[]
gauss_points_results=["X_PHYS"]
VolumeOutput = True
GiDPostMode = "Binary"
GiDWriteMeshFlag = False
GiDWriteConditionsFlag = True
GiDWriteParticlesFlag = False
GiDMultiFileFlag = "Single"
       
# Read optimized model part from restart file
optimized_model_part = ModelPart("optimized_model_part")
optimized_model_part.AddNodalSolutionStepVariable(NORMAL)
restart_file_name = "Small_Cantilever_Restart_File_21"
model_part_io = ModelPartIO(restart_file_name)
model_part_io.ReadModelPart(optimized_model_part)

# Extract volume
threshold = 0.2
extracted_volume_model_part = ModelPart("extracted_volume_model_part")
TopologyExtractorUtilities().ExtractVolumeMesh(optimized_model_part, threshold, extracted_volume_model_part)

# Write extracted volume
gid_io = GiDOutput("extracted_volume_model_part", VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
gid_io.initialize_results(extracted_volume_model_part)
gid_io.write_results(1, extracted_volume_model_part, nodal_results, gauss_points_results)
gid_io.finalize_results()

# Extract surface
extracted_surface_model_part = ModelPart("extracted_surface_model_part")
TopologyExtractorUtilities().ExtractSurfaceMesh(extracted_volume_model_part, extracted_surface_model_part)

# Write extracted surface
gid_io_2 = GiDOutput("extracted_surface_model_part", VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
gid_io_2.initialize_results(extracted_surface_model_part)
gid_io_2.write_results(1, extracted_surface_model_part, nodal_results, gauss_points_results)
gid_io_2.finalize_results()

# Smooth extracted surface
smoothing_relaxation_factor = 1
smoothing_iterations = 5
smoothed_surface_model_part = ModelPart("extracted_surface_model_part")
TopologySmoothingUtilities().SmoothMesh(extracted_surface_model_part, smoothing_relaxation_factor, smoothing_iterations)

# Write smoothed mesh
gid_io_3 = GiDOutput("smoothed_surface_model_part", VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
gid_io_3.initialize_results(extracted_surface_model_part)
gid_io_3.write_results(1, extracted_surface_model_part, nodal_results, gauss_points_results)
gid_io_3.finalize_results()

# Write stl of extracted surface
IOUtilities().WriteSurfaceAsSTLFile("smoothed_design_21.stl",extracted_surface_model_part)