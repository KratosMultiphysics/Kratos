import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

def CreateLineSearch(parameters: Kratos.Parameters):
    type = parameters["type"].GetString()
    if type == "const_step":
        return ConstStep(parameters)
    else:
        raise RuntimeError(f"CreateConvergenceCriteria: unsupported convergence type {type}.")

class ConstStep(object):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"              : "const_step",
            "init_step"          : 0,
        }""")

    def __init__(self, parameters: Kratos.Parameters):
        self.init_step = parameters["init_step"].GetDouble()

    def ComputeStep(self, search_direction):
        norm = KratosOA.ContainerExpressionUtils.NormInf(search_direction)
        if norm:
            step = self.init_step / norm
        else:
            step =  self.init_step
        msg = f"""\t Line Search info: 
            type          : constant 
            value         : {step:0.6e}"""
        print(msg)
        return step