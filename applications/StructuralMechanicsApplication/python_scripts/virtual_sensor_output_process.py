import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils

import os
import pandas as pd
import numpy as np
from numpy import linalg as npla

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VirtualSensorOutputProcess(model, settings["Parameters"])

class VirtualSensorOutputProcess(KratosMultiphysics.OutputProcess):
    def __init__(self, model, settings):
        super().__init__()
        
        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name": "Structure.PLEASE_SPECIFY_SUB_MODEL_PART_NAME",
                "output_control_type": "step",
                "output_interval": 1,
                "intermediate_file_format" : "binary",
                "output_path": "PLEASE_SPECIFY_SUB_MODEL_PART_NAME",
                "save_output_files_in_folder": true,
                "CSV_file_path" : "PLEASE_SPECIFY_CSV_FILE_PATH",
                "CSV_label_template" : "PLEASE_SPECIFY_CSV_LABEL_TEMPLATE",
                "coordinates" : [],
                "nodal_solution_step_data_variables": [
                    "DISPLACEMENT"
                ],
                "nodal_data_value_variables": [],
                "element_data_value_variables": [],
                "condition_data_value_variables": [],
                "gauss_point_variables_extrapolated_to_nodes": []
            }
            """)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[self.settings["model_part_name"].GetString()]

        if self.settings["save_output_files_in_folder"].GetBool():
            if self.model_part.GetCommunicator().MyPID() == 0:
                output_path = settings["output_path"].GetString()
                if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                    kratos_utils.DeleteDirectoryIfExisting(output_path)
                    os.makedirs(output_path)
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        self.next_output = 0.0
        self.__ScheduleNextOutput() # required here esp for restart

        self.map_variable_names_to_variables = {
            "DISPLACEMENT" : KratosMultiphysics.DISPLACEMENT,
            "ACCELERATION" : KratosMultiphysics.ACCELERATION,
            "ROTATION"      : KratosMultiphysics.ROTATION
        }
        
    def PrintOutput(self):

        for nodal_solution_step_data_variable in self.settings["nodal_solution_step_data_variables"]:
            variable_name = nodal_solution_step_data_variable.GetString()
            kratos_variable = self.map_variable_names_to_variables[variable_name]
            
            coordinates_param = self.settings["coordinates"]
            coordinates = []
            for coordinate in coordinates_param:
                coordinates.append(coordinate.GetDouble())

            distances = []
            nodal_values = []
            for node in self.model_part.Nodes:
                node_distance = npla.norm([node.X - coordinates[0], node.Y - coordinates[1], node.Z - coordinates[2]])
                distances.append(node_distance)
                nodal_values.append(node.GetSolutionStepValue(kratos_variable))
            sum_distances = np.sum(distances)
            nodal_weights = np.array(distances)/sum_distances

            node_ctr = 0
            sensor_value = np.zeros(3)
            for node in self.model_part.Nodes:
                sensor_value += nodal_weights[node_ctr] * nodal_values[node_ctr]
                node_ctr += 1
            
            output_path = self.settings["output_path"].GetString()
            sub_model_part_name = self.settings["model_part_name"].GetString().split(".")[-1]
            output_file_path = os.path.join(output_path,sub_model_part_name+"_"+variable_name[0:3]+"_"+str(self.model_part.ProcessInfo[KratosMultiphysics.STEP]))+".npy"
            with open(output_file_path, "wb") as sensor_step_value_file:
                np.save(sensor_step_value_file, sensor_value)

        self.__ScheduleNextOutput()

    def IsOutputStep(self):
        if self.settings["output_control_type"].GetString() == "time":
            return self.__GetTime() >= self.next_output
        else:
            return self.model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.next_output

    def ExecuteFinalize(self):
        super().ExecuteFinalize()

        csv_label_template = self.settings["CSV_label_template"].GetString()

        number_of_steps = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

        num_dimensions = 3
        dimension_labels = ["X","Y","Z"]

        csv_file_path = self.settings["CSV_file_path"].GetString()

        csv_df = pd.read_csv(csv_file_path, delimiter='\t')

        for nodal_solution_step_data_variable in self.settings["nodal_solution_step_data_variables"]:
            variable_name = nodal_solution_step_data_variable.GetString()

            output_path = self.settings["output_path"].GetString()
            sub_model_part_name = self.settings["model_part_name"].GetString().split(".")[-1]
            output_file_path_prefix = os.path.join(output_path,sub_model_part_name+"_"+variable_name[0:3]+"_")
            #str(self.model_part.ProcessInfo[KratosMultiphysics.STEP]))+".npy"
            result_matrix = np.zeros((number_of_steps,num_dimensions))
            for i_step in range(0,number_of_steps):
                output_file_path = output_file_path_prefix + str(i_step+1) + ".npy"
                with open(output_file_path, "rb") as sensor_step_value_file:
                    result_matrix[i_step,:] = np.load(sensor_step_value_file)

            num_missing_values = len(csv_df)-np.shape(result_matrix)[0]
            for i_dim in range(num_dimensions):
                csv_label = variable_name[0:3]+"_"+csv_label_template+"_S_0"+dimension_labels[i_dim]+"0_0"
                csv_df[csv_label] = list(result_matrix[:,i_dim]) + num_missing_values*[None]
        
        new_csv_file_path = csv_file_path.split(".csv")[0]+"_sim.csv"

        csv_df.to_csv(new_csv_file_path, sep='\t', index=False)

        #raise NotImplementedError("This function is not implemented yet!")
        
        #TODO: implement this function to put everything back into the CSV file

    def __ScheduleNextOutput(self):
        if self.settings["output_interval"].GetDouble() > 0.0: # Note: if == 0, we'll just always print
            if self.settings["output_control_type"].GetString() == "time":
                while self.next_output <= self.__GetTime():
                    self.next_output += self.settings["output_interval"].GetDouble()
            else:
                while self.next_output <= self.model_part.ProcessInfo[KratosMultiphysics.STEP]:
                    self.next_output += self.settings["output_interval"].GetDouble()

    def __GetTime(self):
        # remove rounding errors that mess with the comparison
        # e.g. 1.99999999999999999 => 2.0
        return float("{0:.12g}".format(self.model_part.ProcessInfo[KratosMultiphysics.TIME]))