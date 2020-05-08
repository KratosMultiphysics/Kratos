from KratosMultiphysics import *
import KratosMultiphysics.ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

from math import *

class TestSymbolicEulerianConvectionDiffusionElement(ConvectionDiffusionAnalysis):
    def __init__(self,model,parameters):
        super(TestSymbolicEulerianConvectionDiffusionElement,self).__init__(model,parameters)
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
            x = node.X
            y = node.Y
            diffusivity = 2.0 # check in materials.json
            # forcing = -432.0 * (x**2 + y**2 - x - y)
            forcing = ( 4*pi*x*(-sin(2*pi*x)*sin(2*pi*y)+cos(2*pi*x)) + 6*pi*y*cos(2*pi*x)*cos(2*pi*y) ) + \
                diffusivity*( -4*(pi)**2*cos(2*pi*x)*sin(2*pi*y) - 4*(pi)**2*sin(2*pi*x) - 4*(pi)**2*cos(2*pi*x)*sin(2*pi*y) )
            convective_velocity = [2.0*x,3.0*y,0.0]
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,convective_velocity)
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing)
        for node in self.model.GetModelPart("MainModelPart.Subpart_Boundary").Nodes:
            x = node.X
            y = node.Y
            temperature = cos(2*pi*x)*sin(2*pi*y)#+sin(2*pi*x)
            print(x,y,temperature)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,temperature)



if __name__ == "__main__":
    from sys import argv

    project_parameters_file_name = "test_symbolic_eulerian_convection_diffusion_explicit_element/project_parameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = TestSymbolicEulerianConvectionDiffusionElement(model, parameters)
    simulation.Run()