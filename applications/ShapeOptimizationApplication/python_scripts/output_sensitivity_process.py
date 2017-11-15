import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as ShapeOptimizationApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return OutputSensitivityProcess(Model, settings["Parameters"])

class OutputSensitivityProcess(ShapeOptimizationApplication.OutputSensitivityProcess):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]  
        
        ShapeOptimizationApplication.OutputSensitivityProcess.__init__(self, model_part, settings)
