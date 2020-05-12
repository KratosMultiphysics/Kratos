from KratosMultiphysics import *
import KratosMultiphysics.ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

from math import *

class TestSymbolicEulerianConvectionDiffusionElement(ConvectionDiffusionAnalysis):
    def __init__(self,model,parameters):
        super(TestSymbolicEulerianConvectionDiffusionElement,self).__init__(model,parameters)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)

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
            x = node.X
            y = node.Y
            diffusivity = 2.0 # check in materials.json
            convective_velocity = [2.0,3.0,0.0]

            forcing = -diffusivity*(-4*pi*pi*cos(2*pi*x)*sin(2*pi*y)-4*pi*pi*cos(2*pi*x)*sin(2*pi*y)) + \
                convective_velocity[0]*(-2*pi*sin(2*pi*x)*sin(2*pi*y)) + convective_velocity[1]*(2*pi*cos(2*pi*x)*cos(2*pi*y))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,convective_velocity)
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing)

if __name__ == "__main__":
    from sys import argv

    project_parameters_file_name = "test_symbolic_eulerian_convection_diffusion_explicit_element/project_parameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = TestSymbolicEulerianConvectionDiffusionElement(model, parameters)
    simulation.Run()

    # check L2 error via midpoint rule and estimate minimal nodal h
    KratosMultiphysics.CalculateNodalAreaProcess(simulation._GetSolver().main_model_part,2).Execute()
    find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(simulation._GetSolver().main_model_part)
    find_nodal_h.Execute()
    error = 0
    nodal_h = 1e9
    model_part_name = simulation.project_parameters["problem_data"]["model_part_name"].GetString()
    for node in simulation.model.GetModelPart(model_part_name).Nodes:
        # L2 norm
        x = node.X
        y = node.Y
        u_analytical = cos(2*pi*x)*sin(2*pi*y)
        u_numerical = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
        error = error + ((u_analytical - u_numerical)**2*node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        # nodal h
        if (node.GetValue(KratosMultiphysics.NODAL_H) < nodal_h):
            nodal_h = node.GetValue(KratosMultiphysics.NODAL_H)

    error = sqrt(error)
    print("minimal nodal h:",nodal_h)
    print("error:",error)
