import KratosMultiphysics
import KratosMultiphysics.PfemApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")

    return ConstantRotationProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ConstantRotationProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)


    def ExecuteInitialize(self):
        name_of_submodelpart = settings["model_part_name"].GetString()
        self.model_part = Model[name_of_submodelpart]
        self.cplusplus_version_process = KratosMultiphysics.PfemApplication.ConstantRotationProcess(self.model_part, settings)


    def ExecuteInitializeSolutionStep(self):
        self.cplusplus_version_process.ExecuteInitializeSolutionStep()
