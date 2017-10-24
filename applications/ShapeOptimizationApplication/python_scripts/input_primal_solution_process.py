import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as ShapeOptimizationApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return InputPrimalSolutionProcess(Model, settings["Parameters"])

class InputPrimalSolutionProcess(ShapeOptimizationApplication.InputPrimalSolutionProcess):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        
        ShapeOptimizationApplication.InputPrimalSolutionProcess.__init__(self, self.model_part, settings)

# fusseder TODO: file copied, maybe things have to be renamed