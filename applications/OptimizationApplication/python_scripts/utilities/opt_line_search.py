import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def CreateLineSearch(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    type = parameters["type"].GetString()
    if type == "const_step":
        return ConstStep(parameters, optimization_problem)
    else:
        raise RuntimeError(f"CreateConvergenceCriteria: unsupported convergence type {type}.")

class ConstStep(object):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"              : "const_step",
            "init_step"          : 0,
            "gradient_scaling": "inf_norm"
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.init_step = parameters["init_step"].GetDouble()
        self.__optimization_problem = optimization_problem

    def ComputeStep(self) -> float:
        algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()
        
        if not algorithm_buffered_data.HasValue("search_direction"):
            raise RuntimeError(f"Algorithm data does not contain computed \"search_direction\".\nData:\n{algorithm_buffered_data}" )
        
        norm = KratosOA.ContainerExpressionUtils.NormInf(algorithm_buffered_data["search_direction"])
        if norm:
            step = self.init_step / norm
        else:
            step =  self.init_step
        msg = f"""\t Line Search info: 
            type          : constant 
            value         : {step:0.6e}"""
        print(msg)
        return step