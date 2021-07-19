import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
import numpy as np
import os

def Factory(parameters, model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return RomBasisProcess(model, parameters["Parameters"])

## All the processes python should be derived from "Process"

class RomBasisProcess(KratosMultiphysics.Process):
    def __init__(self, Model, parameters):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "main_model_part": "Structure",
            "time_steps_per_file": 100,
            "folder_name": "SnapshotsMatrix"
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        parameters.ValidateAndAssignDefaults(default_settings)
        # Alternative:
        # settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.parameters = parameters
        self.model_part_name = self.parameters["main_model_part"].GetString()
        self.model_part = Model[self.model_part_name]
        self.snapshots_matrix = []
        self.time_steps_per_file = self.parameters["time_steps_per_file"].GetInt()
        self.folder_name = self.parameters["folder_name"].GetString()
        self.parameters.RemoveValue("time_steps_per_file")
        self.parameters.RemoveValue("folder_name")
        
        
    def ExecuteInitialize(self):
        self.output_utility = romapp.RomOutputUtility(self.model_part,self.parameters)
        self.output_utility.Start()
        self.counter = 0
        self.name = 0
        try:
            os.mkdir(self.folder_name)
        except:
            pass
    
    def ExecuteFinalizeSolutionStep(self):
        if self.model_part_name=="Structure":
            self.snapshots_matrix.append(self.output_utility.GetSnapshotOfDisplacements())
            self.counter += 1
        elif self.model_part_name=="FluidModelPart":
            self.snapshots_matrix.append(self.output_utility.GetSnapshotOfVelocityAndPressure())
            self.counter += 1
        else:
            raise Exception("model_part_name should be Structure or FluidModelPart")
        if self.counter==self.time_steps_per_file:
            snapshots_matrix_array = np.array(self.snapshots_matrix,copy=False).T 
            self.name += 1
            with open(os.path.join(self.folder_name,'SnapshotsMatrix_'+str(self.name)+'.npy'), 'wb') as f:
                np.save(f, snapshots_matrix_array)
                self.snapshots_matrix = []
                self.counter = 0
    
    def ExecuteFinalize(self):
        if self.counter==0:
            pass
        else:
            snapshots_matrix_array = np.array(self.snapshots_matrix)
            self.name += 1
            with open(os.path.join(self.folder_name,'SnapshotsMatrix_'+str(self.name)+'.npy'), 'wb') as f:
                np.save(f, snapshots_matrix_array)
