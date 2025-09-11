from kratos_external_solver_optimization import external_optimization
from dummy_solver import DummySolver

external_optimization = external_optimization.ExternalOptimization(
    "optimization_parameters.json",
    DummySolver,
    "dummy_solver_parameters.json")

external_optimization.RunOptimization()
