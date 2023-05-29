import KratosMultiphysics as KM
import numpy as np
import os

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VariableSnapshotProcess(model, settings["Parameters"])

class VariableSnapshotProcess(KM.Process):

    def __init__(self, model, settings):
        KM.Process.__init__(self)  # calling the baseclass constructor

        # define default settings
        default_settings = KM.Parameters("""{
            "model_part_name" : "please_specify_model_part_name",
            "list_of_variables" : [],
            "file_prefix" : "",
            "folder_name" : ""
        }""")
        
        # validate and assign default settings
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.num_nodes = len(self.model_part.Nodes)
        self.variables = settings["list_of_variables"].GetStringArray()
        self.prefix = settings["file_prefix"].GetString()
        self.folder_name = settings["folder_name"].GetString()

        # If folder_name is not an empty string, create the folder
        if self.folder_name:
            os.makedirs(self.folder_name, exist_ok=True)

        # Initialize a dictionary to hold data for each variable
        self.data = {var: [] for var in self.variables}

        # Create mapping from node ID to index
        self.node_id_to_index = {node.Id: i for i, node in enumerate(self.model_part.Nodes)}

    def ExecuteFinalizeSolutionStep(self):
        for var_name in self.variables:
            var = KM.KratosGlobals.GetVariable(var_name)  # Get the variable
            data_step = np.zeros(self.num_nodes)  # Initialize array for this step
            
            for node_id, i in self.node_id_to_index.items():
                node = self.model_part.Nodes[node_id]
                data_step[i] = node.GetSolutionStepValue(var)
            
            self.data[var_name].append(data_step)

    def ExecuteFinalize(self):
        for var_name in self.variables:
            data_all_steps = np.stack(self.data[var_name])

            # If folder_name is not an empty string, prepend it to the file name
            if self.folder_name:
                file_name = f"{self.folder_name}/{self.prefix}_{var_name}.npy"
            else:
                file_name = f"{self.prefix}_{var_name}.npy"

            np.save(file_name, data_all_steps)

        # Save the node ID to index mapping as well
        if self.folder_name:
            file_name = f"{self.folder_name}/{self.prefix}_node_id_to_index.npy"
        else:
            file_name = f"{self.prefix}_node_id_to_index.npy"

        np.save(file_name, self.node_id_to_index)

