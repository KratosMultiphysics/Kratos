import KratosMultiphysics as km
from KratosMultiphysics.TopologyOptimizationApplication import topology_optimizer_factory

with open("TopologyOptimizationParameters.json",'r') as parameter_file:
    parameters = km.Parameters(parameter_file.read())
model = km.Model()

topology_optimizer = topology_optimizer_factory.CreateTopologyOptimizer(parameters["topology_optimization_settings"], model)
topology_optimizer.Optimize()
topology_optimizer.DeleteVoidElements()
topology_optimizer.Reanalyze()
topology_optimizer.Smooth()
topology_optimizer.Reanalyze()
topology_optimizer.ExportStl()
topology_optimizer.ExportStp()
