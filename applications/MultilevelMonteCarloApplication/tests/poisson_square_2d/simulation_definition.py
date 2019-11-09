from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.MultilevelMonteCarloApplication as KratosMLMC
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

# Importing the problem analysis stage class
from  KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

# Avoid printing of Kratos informations
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

"""
SimulationScenario is inherited from the Analysis Stage class and solves the Poisson PDE in domain \Omega = (0,1)^2 with zero Dirichlet boundary conditions
-lapl(u) = \varepsilon*f    u \in \Omega
u = 0                       u \in \partial(\Omega)
where
\varepsilon \sim Beta(2,6)
and f = f1 or f = f2, with
f1= -432*x*(x-1)*y*(y-1)
f2= -432*(x**2+y**2-x-y)
the QoI is
QoI = \int_(\Omega) u(x,y)dxdy
"""
class SimulationScenario(ConvectionDiffusionAnalysis):
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        super(SimulationScenario,self).__init__(input_model,input_parameters)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _CreateSolver(self):
        from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_stationary_solver
        return convection_diffusion_stationary_solver.CreateSolver(self.model,self.project_parameters["solver_settings"])

    """
    function introducing the stochasticity in the right hand side
    input:  self: an instance of the class
    """
    def ModifyInitialProperties(self):
        model_part_name = self.project_parameters["problem_data"]["model_part_name"].GetString()
        for node in self.model.GetModelPart(model_part_name).Nodes:
            coord_x = node.X
            coord_y = node.Y
            # forcing = -432.0 * coord_x * (coord_x - 1) * coord_y * (coord_y - 1)
            forcing = -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y)
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing*self.sample[0])

    """
    function evaluating the QoI of the problem: int_{domain} TEMPERATURE(x,y) dx dy
    midpoint rule used to compute the integral
    input:  self: an instance of the class
    """
    def EvaluateQuantityOfInterest(self):
        KratosMultiphysics.CalculateNodalAreaProcess(self._GetSolver().main_model_part,2).Execute()
        Q = 0.0
        for node in self._GetSolver().main_model_part.Nodes:
            Q = Q + (node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)*node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        return Q