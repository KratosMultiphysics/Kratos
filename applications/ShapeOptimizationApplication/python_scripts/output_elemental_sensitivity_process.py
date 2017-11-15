import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as ShapeOptimizationApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return OutputElementalSensitivityProcess(Model, settings["Parameters"])

class OutputElementalSensitivityProcess(ShapeOptimizationApplication.OutputElementalSensitivityProcess):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]  
        
        ShapeOptimizationApplication.OutputElementalSensitivityProcess.__init__(self, model_part, settings)
