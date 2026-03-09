import KratosMultiphysics
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ControlModuleProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ControlModuleProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        # Control module process acting on the imposed direction: 0 (X), 1 (Y), 2 (Z) or 3 (radial)
        # The radial direction is valid only for the vertical walls of a right cylinder with the base
        # on the 'X-Y' plane centered on (0,0). Negative target_stress means compression.

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.control_module_process = DemFem.ControlModuleProcess(self.model_part, settings)

    def ExecuteInitialize(self):
        self.control_module_process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.control_module_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.control_module_process.ExecuteFinalizeSolutionStep()