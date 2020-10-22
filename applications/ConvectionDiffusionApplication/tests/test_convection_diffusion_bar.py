from KratosMultiphysics import *
import KratosMultiphysics.ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import *

class TestConvectionDiffusionBar(KratosUnittest.TestCase):

    def testConvectionDiffusionBarExplicitElementUnsteadyDOSS(self):
        project_parameters_file_name = "test_convection_diffusion_bar/project_parameters_bar_DOSS.json"
        with open(project_parameters_file_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        model = KratosMultiphysics.Model()
        bar_simulation = ConvectionDiffusionAnalysis(model, parameters)
        bar_simulation._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        bar_simulation.Run()
        # check L2 error via midpoint rule
        KratosMultiphysics.CalculateNodalAreaProcess(bar_simulation._GetSolver().main_model_part,2).Execute()
        error = 0
        model_part_name = bar_simulation.project_parameters["problem_data"]["model_part_name"].GetString()
        for node in bar_simulation.model.GetModelPart(model_part_name).Nodes:
            # L2 norm
            x = node.X
            u_analytical = x - sin(bar_simulation.time)
            u_numerical = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            error += (((u_analytical - u_numerical)**2)*node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        error = sqrt(error)
        self.assertAlmostEqual(error,0.0017678016487118946,delta=1e-12)

    def testConvectionDiffusionBarExplicitElementUnsteadyQOSS(self):
        project_parameters_file_name = "test_convection_diffusion_bar/project_parameters_bar_QOSS.json"
        with open(project_parameters_file_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        model = KratosMultiphysics.Model()
        bar_simulation = ConvectionDiffusionAnalysis(model, parameters)
        bar_simulation._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        bar_simulation.Run()
        # check L2 error via midpoint rule
        KratosMultiphysics.CalculateNodalAreaProcess(bar_simulation._GetSolver().main_model_part,2).Execute()
        error = 0
        model_part_name = bar_simulation.project_parameters["problem_data"]["model_part_name"].GetString()
        for node in bar_simulation.model.GetModelPart(model_part_name).Nodes:
            # L2 norm
            x = node.X
            u_analytical = x - sin(bar_simulation.time)
            u_numerical = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            error += (((u_analytical - u_numerical)**2)*node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        error = sqrt(error)
        self.assertAlmostEqual(error,0.0017678016487242358,delta=1e-12)

    def testConvectionDiffusionBarExplicitElementUnsteadyDASGS(self):
        project_parameters_file_name = "test_convection_diffusion_bar/project_parameters_bar_DASGS.json"
        with open(project_parameters_file_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        model = KratosMultiphysics.Model()
        bar_simulation = ConvectionDiffusionAnalysis(model, parameters)
        bar_simulation._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        bar_simulation.Run()
        # check L2 error via midpoint rule
        KratosMultiphysics.CalculateNodalAreaProcess(bar_simulation._GetSolver().main_model_part,2).Execute()
        error = 0
        model_part_name = bar_simulation.project_parameters["problem_data"]["model_part_name"].GetString()
        for node in bar_simulation.model.GetModelPart(model_part_name).Nodes:
            # L2 norm
            x = node.X
            u_analytical = x - sin(bar_simulation.time)
            u_numerical = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            error += (((u_analytical - u_numerical)**2)*node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        error = sqrt(error)
        self.assertAlmostEqual(error,0.04953302968545053,delta=1e-12)

    def testConvectionDiffusionBarExplicitElementUnsteadyQASGS(self):
        project_parameters_file_name = "test_convection_diffusion_bar/project_parameters_bar_QASGS.json"
        with open(project_parameters_file_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        model = KratosMultiphysics.Model()
        bar_simulation = ConvectionDiffusionAnalysis(model, parameters)
        bar_simulation._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        bar_simulation.Run()
        # check L2 error via midpoint rule
        KratosMultiphysics.CalculateNodalAreaProcess(bar_simulation._GetSolver().main_model_part,2).Execute()
        error = 0
        model_part_name = bar_simulation.project_parameters["problem_data"]["model_part_name"].GetString()
        for node in bar_simulation.model.GetModelPart(model_part_name).Nodes:
            # L2 norm
            x = node.X
            u_analytical = x - sin(bar_simulation.time)
            u_numerical = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            error += (((u_analytical - u_numerical)**2)*node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        error = sqrt(error)
        self.assertAlmostEqual(error,0.045232415452161175,delta=1e-12)

if __name__ == '__main__':
    KratosUnittest.main()

