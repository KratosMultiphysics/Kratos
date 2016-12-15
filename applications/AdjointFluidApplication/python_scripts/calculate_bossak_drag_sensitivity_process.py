import KratosMultiphysics
import KratosMultiphysics.AdjointFluidApplication as AdjointFluidApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CalculateBossakDragSensitivityProcess(Model, settings["Parameters"])

class CalculateBossakDragSensitivityProcess(AdjointFluidApplication.CalculateBossakDragSensitivityProcess):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]  
        
        AdjointFluidApplication.CalculateBossakDragSensitivityProcess.__init__(self, model_part, settings)
