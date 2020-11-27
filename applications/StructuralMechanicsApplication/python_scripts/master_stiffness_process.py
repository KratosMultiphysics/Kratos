import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def Factory(parameters, model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MasterStiffnessProcess(model, parameters["Parameters"])

## All the processes python should be derived from "Process"

class MasterStiffnessProcess(KratosMultiphysics.Process):
    def __init__(self, Model, parameters):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model["Structure"]
        self.master_stiffness_process = StructuralMechanicsApplication.MasterStiffnessProcess(self.model_part,parameters)

    # This method is executed in order to initialize the current step
    def ExecuteInitializeSolutionStep(self):
        self.master_stiffness_process.ExecuteInitializeSolutionStep()

    # This method is executed in order to finalize the current step
    def ExecuteFinalizeSolutionStep(self):
        self.master_stiffness_process.ExecuteFinalizeSolutionStep()