# Import Kratos core and apps
import os, shutil
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import structural_response_function_factory
import structural_mechanics_analysis
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosHDF5Application
import h5py
import numpy as np
import matplotlib.pyplot as plt
import time as timer

with open("nonlinear_shell_test_parameters.json",'r') as parameter_file:
    parameters = KratosMultiphysics.Parameters( parameter_file.read())

## input values
number_of_nodes = 4
number_of_dofs = 24

model = KratosMultiphysics.Model()

model_part_primal = model.CreateModelPart("Structure" , 2)
primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, parameters)
primal_analysis.Initialize()

startTime = timer.time()
if not primal_analysis.time < primal_analysis.end_time:
    primal_analysis.end_time += 1

primal_analysis.RunSolutionLoop()

file_name = []
for time_step in np.arange(0.1, 5.01, 0.1):
    time_step_string = f"{time_step:.4f}"
    file_name.append("primal_output_shell-" + time_step_string + ".h5") 
    
element_connectivities = []
for element in model_part_primal.Elements:
    Nodes = element.GetNodes()
    node_Ids = []
    for node in Nodes:
        node_Ids.append(node.Id)
    element_connectivities.append(node_Ids)

## reading the hdf files data
elements_displacement_vectors = []
elements_load_vectors = []

## concatanate the an array at time t = 0
elements_displacement_vectors.append(np.zeros(number_of_dofs))
elements_load_vectors.append(np.zeros(number_of_dofs))

number_of_time_steps = 49
for i in range(0,number_of_time_steps+1):
    hdf5_file = h5py.File(file_name[i], 'r')
    nodal_data = hdf5_file["ResultsData/NodalSolutionStepData"]
    nodal_displacement_hdf = np.array(nodal_data['DISPLACEMENT'])
    nodal_point_load_hdf = np.array(nodal_data['POINT_LOAD'])
    nodal_reaction_force_hdf = np.array(nodal_data['REACTION'])
    hdf5_file.close()
    element_displacement_single_step = []
    element_load_single_step = []
    for element in model_part_primal.Elements:
        node_displacement = []
        node_load = []
        for i in range(0,number_of_nodes):
            node_displacement.append(nodal_displacement_hdf[element_connectivities[element.Id - 1][i] - 1])
            node_load.append(nodal_point_load_hdf[element_connectivities[element.Id - 1][i] - 1])
        ## shell elements needs to concatanate values at the four nodes
        element_displacement = np.concatenate([node_displacement[0], node_displacement[1] , node_displacement[2], node_displacement[3] ])
        element_displacement_single_step.append(element_displacement)
        element_load = np.concatenate([node_load[0], node_load[1] , node_load[2], node_load[3] ])
        element_load_single_step.append(element_load)    
    
    elements_displacement_vectors.append(element_displacement_single_step)
    elements_load_vectors.append(element_load_single_step)
print("displacement vector: " , elements_displacement_vectors[49])
print("laod vector: " , elements_load_vectors[49])

## calculating the elastic strain energy
response_value = 0
for element in model_part_primal.Elements:
    element_response = 0
    for i in range(1, number_of_time_steps + 2):
        disp_step = elements_displacement_vectors[i][element.Id - 1] - elements_displacement_vectors[i - 1][element.Id - 1]
        load = 0.5 * (elements_load_vectors[i][element.Id - 1] + elements_load_vectors[i - 1][element.Id - 1]) / 2
        element_response = element_response + np.inner(disp_step , load)
    print("element response",element_response)
    response_value = response_value + element_response
print("strain energy" , response_value)


###### reading and outputting the data to kratos
# file_name = []
# for time_step in np.arange(0.1 ,1.01, 0.1):
#     time_step_string = f"{time_step:.4f}"
#     file_name.append("primal_output_truss-" + time_step_string + ".h5") 

# nodal_displacement_x = []
# nodal_displacement_y = []
# nodal_displacement_z = []
# nodal_point_load_x = []
# nodal_point_load_y = []
# nodal_point_load_z = []
# nodal_reaction_force_x = []
# nodal_reaction_force_y = []
# nodal_reaction_force_z = []

# number_of_time_steps = 9

# ## reading the hdf files data
# for i in range(0,number_of_time_steps+1):
#     hdf5_file = h5py.File(file_name[i], 'r')
#     nodal_data = hdf5_file["ResultsData/NodalSolutionStepData"]
#     nodal_displacement_hdf = np.array(nodal_data['DISPLACEMENT'])
#     nodal_point_load_hdf = np.array(nodal_data['POINT_LOAD'])
#     nodal_reaction_force_hdf = np.array(nodal_data['REACTION'])
#     nodal_displacement_x.append(nodal_displacement_hdf[:,0].transpose())
#     nodal_displacement_y.append(nodal_displacement_hdf[:,1].transpose())
#     nodal_displacement_z.append(nodal_displacement_hdf[:,2].transpose())
#     nodal_point_load_x.append(nodal_point_load_hdf[:,0].transpose())
#     nodal_point_load_y.append(nodal_point_load_hdf[:,1].transpose())
#     nodal_point_load_z.append(nodal_point_load_hdf[:,2].transpose())
#     nodal_reaction_force_x.append(nodal_reaction_force_hdf[:,0].transpose())
#     nodal_reaction_force_y.append(nodal_reaction_force_hdf[:,1].transpose())
#     nodal_reaction_force_z.append(nodal_reaction_force_hdf[:,2].transpose())
#     hdf5_file.close()

# write out another hdf5 file with the data over all the time steps
# shutil.copyfile(file_name[0], "primal_output_truss_total-1.0000.h5")
# hdf_output = h5py.File("primal_output_truss_total-1.0000.h5", "a")
# hdf_output.create_group("ResultsData/NodalDataValues")
# hdf_output["ResultsData/NodalDataValues/DISPLACEMENT_X_TOTAL"] = np.array(nodal_displacement_x).transpose()
# hdf_output["ResultsData/NodalDataValues/DISPLACEMENT_Y_TOTAL"] = np.array(nodal_displacement_y).transpose()
# hdf_output["ResultsData/NodalDataValues/DISPLACEMENT_Z_TOTAL"] = np.array(nodal_displacement_z).transpose()
# hdf_output["ResultsData/NodalDataValues/POINT_LOAD_X_TOTAL"] = np.array(nodal_point_load_x).transpose()
# hdf_output["ResultsData/NodalDataValues/POINT_LOAD_Y_TOTAL"] = np.array(nodal_point_load_y).transpose()
# hdf_output["ResultsData/NodalDataValues/POINT_LOAD_Z_TOTAL"] = np.array(nodal_point_load_z).transpose()
# hdf_output["ResultsData/NodalDataValues/REACTION_X_TOTAL"] = np.array(nodal_reaction_force_x).transpose()
# hdf_output["ResultsData/NodalDataValues/REACTION_Y_TOTAL"] = np.array(nodal_reaction_force_y).transpose()
# hdf_output["ResultsData/NodalDataValues/REACTION_Z_TOTAL"] = np.array(nodal_reaction_force_z).transpose()
# hdf_output["ResultsData/NodalDataValues/DISPLACEMENT_X_TOTAL"].attrs.create("Size1", number_of_time_steps+1)
# hdf_output["ResultsData/NodalDataValues/DISPLACEMENT_Y_TOTAL"].attrs.create("Size1", number_of_time_steps+1)
# hdf_output["ResultsData/NodalDataValues/DISPLACEMENT_Z_TOTAL"].attrs.create("Size1", number_of_time_steps+1)
# hdf_output["ResultsData/NodalDataValues_partition"] = list(hdf_output["ResultsData/NodalSolutionStepData_partition"])
# hdf_output.close()



# output the results to GID
# from gid_output_process import GiDOutputProcess
# gid_output = GiDOutputProcess(model_part_primal,
#                              "gid_output",
#                             KratosMultiphysics.Parameters("""
#                                 {
#                                     "result_file_configuration" : {
#                                         "gidpost_flags": {
#                                             "GiDPostMode": "GiD_PostBinary",
#                                             "WriteDeformedMeshFlag": "WriteUndeformed",
#                                             "WriteConditionsFlag": "WriteConditions",
#                                             "MultiFileFlag": "SingleFile"
#                                         },
#                                         "nodal_results"       : ["DISPLACEMENT"],
#                                         "output_frequency"    : 1
#                                     }
#                                 }
#                                 """)
#                             )

# gid_output.ExecuteInitialize()
# gid_output.ExecuteBeforeSolutionLoop()
# gid_output.ExecuteInitializeSolutionStep()
# gid_output.PrintOutput()
# gid_output.ExecuteFinalizeSolutionStep()
# gid_output.ExecuteFinalize()



