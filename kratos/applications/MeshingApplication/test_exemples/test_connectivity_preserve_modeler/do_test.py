from KratosMultiphysics import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *

#importing the benchmarking auxiliary functions
kratos_benchmarking_path = '../../../../benchmarking'
import sys
sys.path.append(kratos_benchmarking_path)
import benchmarking

#defining a model part
origin_model_part = ModelPart("FluidPart");  
new_model_part = ModelPart("NewPart");

#define variables for the solvers to be used
import fractional_step_solver
fractional_step_solver.AddVariables(origin_model_part)
origin_model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
origin_model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
origin_model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)

##ATTENTION: IMPORTANT!! note that the convection_diffusion_solver will work on a different model part
##nevertheless, since the new model part will be filled up later, and the nodes to be created will 
##be shared between the origin and new model part, we will add the variables to the "origin_model_part"
##such variables will be copied when creating then new one
thermal_settings = ConvectionDiffusionSettings()
thermal_settings.SetDensityVariable(DENSITY)
thermal_settings.SetDiffusionVariable(CONDUCTIVITY)
thermal_settings.SetUnknownVariable(TEMPERATURE)
thermal_settings.SetVolumeSourceVariable(HEAT_FLUX)
thermal_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
thermal_settings.SetMeshVelocityVariable(MESH_VELOCITY)
thermal_settings.SetProjectionVariable(TEMP_CONV_PROJ);

import convection_diffusion_solver
convection_diffusion_solver.AddVariables(origin_model_part,thermal_settings)

##read the original_model part
origin_file_name = "../Mapping_2d.gid/Mapping_2d"
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(origin_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_pfem = ModelPartIO(origin_file_name)
model_part_io_pfem.ReadModelPart(origin_model_part)

buffer_size = 3
origin_model_part.SetBufferSize(buffer_size)

##generate a new property in the origin_model_part to insert control values
prop_initial_value = 123.0
prop_test_value = prop_initial_value
origin_model_part.Properties[10].SetValue(PRESSURE,prop_initial_value)

##generate a new model part with different elements. Will use THE SAME Nodes, Properties and ProcessInfo as the original one
##but will have different Elements and Conditions
connectivity_model_part_generator = ConnectivityPreserveModeler()
connectivity_model_part_generator.GenerateModelPart(origin_model_part,new_model_part,"ConvDiff2D","ThermalFace2D")

#verify buffer size is assigned correctly
benchmarking.Output(origin_model_part.GetBufferSize(), "test expected value for buffer_size on origin_model_part", buffer_size)
benchmarking.Output(new_model_part.GetBufferSize(), "test expected value for buffer_size on new_model_part", buffer_size)

##loop over time and verify that the ProcessInfo behaves correctly
dt = 0.1
time = 0.0
for i in range(0,30):
  time = time + dt
  
  #we advance in time ONE of the model_parts
  #note that it would be an error to advance both of them since the processinfo is shared
  origin_model_part.CloneTimeStep(time)
  
  if(i > 2):
    #here we verify that both model parts have the correct time and delta time
    benchmarking.Output(origin_model_part.ProcessInfo[TIME], "expected time", time)
    benchmarking.Output(origin_model_part.ProcessInfo[DELTA_TIME], "expected delta time", dt)
    benchmarking.Output(new_model_part.ProcessInfo[TIME], "expected time", time)
    benchmarking.Output(new_model_part.ProcessInfo[DELTA_TIME], "expected delta time", dt)
    
    #verify that values of the processinfo can be set correctly in both model parts
    test_value = 1234.0
    origin_model_part.ProcessInfo.SetValue(YOUNG_MODULUS,test_value)    #here we set from origin model part
    benchmarking.Output(origin_model_part.ProcessInfo[YOUNG_MODULUS], "test expected value for YOUNG_MODULUS on origin_model_part", test_value)
    benchmarking.Output(new_model_part.ProcessInfo[YOUNG_MODULUS], "test expected value for YOUNG_MODULUS on new_model_part", test_value)
    
    
    other_test_value = 4567.0
    new_model_part.ProcessInfo.SetValue(TEMPERATURE,other_test_value)    #here we set from new_model_part
    benchmarking.Output(origin_model_part.ProcessInfo[TEMPERATURE], "test expected value for TEMPERATURE  on origin_model_part", other_test_value)
    benchmarking.Output(new_model_part.ProcessInfo[TEMPERATURE], "test expected value for TEMPERATURE  on new_model_part", other_test_value)
    
    ##add 1.0 to the value of PRESSURE for the model part
    ##and verify that this is correct on both model parts
    prop_test_value += 2.0;
    aaa = origin_model_part.Properties[10].GetValue(PRESSURE) 
    origin_model_part.Properties[10].SetValue(PRESSURE,aaa+1.0)    
    aaa = new_model_part.Properties[10].GetValue(PRESSURE) 
    new_model_part.Properties[10].SetValue(PRESSURE,aaa+1.0)    
    benchmarking.Output(origin_model_part.Properties[10].GetValue(PRESSURE) , "test expected value for PRESSURE in properties for the origin_model_part", prop_test_value)
    benchmarking.Output(new_model_part.Properties[10].GetValue(PRESSURE)    , "test expected value for PRESSURE in properties for the new_model_part", prop_test_value)

print "completed"
  
