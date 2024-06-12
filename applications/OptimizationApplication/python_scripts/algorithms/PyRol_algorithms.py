import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return PYROLAlgorithms(model, parameters, optimization_problem)

class PYROLAlgorithms():
    pass