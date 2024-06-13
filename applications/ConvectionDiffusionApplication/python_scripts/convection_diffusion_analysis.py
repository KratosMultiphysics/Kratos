
# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion as solver_wrapper
import KratosMultiphysics.ConstitutiveLawsApplication as CLA

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

# Other imports
import sys

class ConvectionDiffusionAnalysis(AnalysisStage):
    """
    This class is the main-script of the ConvectionDiffusionApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        solver_settings = project_parameters["solver_settings"]

        if not solver_settings.Has("domain_size"):
            KratosMultiphysics.Logger.PrintInfo("ConvectionDiffusionAnalysis", "Using the old way to pass the domain_size, this will be removed!")
            solver_settings.AddEmptyValue("domain_size")
            solver_settings["domain_size"].SetInt(project_parameters["problem_data"]["domain_size"].GetInt())

        super(ConvectionDiffusionAnalysis, self).__init__(model, project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[Convection-Diffusion Simulation]:: "

    def RunSolutionLoop(self):
            """This function executes the solution loop of the AnalysisStage
            It can be overridden by derived classes
            """
            while self.KeepAdvancingSolutionLoop():
                self.time = self._GetSolver().AdvanceInTime(self.time)
                process = CLA.AdvanceInTimeHighCycleFatigueProcess(self._GetSolver().GetComputingModelPart(), self.project_parameters)
                if self.project_parameters.Has("fatigue"):
                    if self.project_parameters["fatigue"]["advancing_strategy"].GetBool():
                        process.Execute()
                        time_incr = self._GetSolver().GetComputingModelPart().ProcessInfo[CLA.TIME_INCREMENT]
                        self.time += time_incr
                        self._GetSolver().GetComputingModelPart().ProcessInfo[CLA.TIME_INCREMENT] = 0.0
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time
                self.InitializeSolutionStep()
                self._GetSolver().Predict()
                self._GetSolver().SolveSolutionStep()
                self.FinalizeSolutionStep()
                self.OutputSolutionStep()

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 convection_diffusion_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 convection_diffusion_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ConvectionDiffusionAnalysis(model, parameters)
    simulation.Run()
