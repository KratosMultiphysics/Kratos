from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
# Importing the base class
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class StructuralMechanicsPrebucklingAnalysis(StructuralMechanicsAnalysis):
    def __init__(self, model, project_parameters):
        super(StructuralMechanicsPrebucklingAnalysis, self).__init__(model, project_parameters)

    def Initialize(self):
        """This function initializes the StructuralMechanicsPrebucklingAnalysis
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        problem_data = self.project_parameters["problem_data"]
        if problem_data.Has("start_time"):
            warn_msg = 'Parameter TIME is used as load factor. \n'
            warn_msg += 'Parameter "start_time" will be ignored!'
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsPrebucklingAnalysis; Warning", warn_msg)
        else:
            # Create dummy parameter
            aux_settings = KratosMultiphysics.Parameters(r"""{ "start_time" : 1.0 }""")
            problem_data.AddMissingParameters(aux_settings)

        if problem_data.Has("end_time"):
            warn_msg = 'Parameter TIME is used as load factor. \n'
            warn_msg += 'Parameter "end_time" will be ignored!'
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsPrebucklingAnalysis; Warning", warn_msg)
        else:
            # Create dummy paramter
            aux_settings = KratosMultiphysics.Parameters(r"""{ "end_time" : 1.0 }""")
            problem_data.AddMissingParameters(aux_settings)

        # Initialize super class
        super().Initialize()

        # Initialize solution stepping
        self.step = 0
        self.time = 1
        if not problem_data.Has("nsteps"):
            raise Exception("StructuralMechanicsPrebucklingAnalysis: " + 'Maximum number of steps "nsteps" must be provided"!')
        else:
            self.nsteps = problem_data["nsteps"].GetInt()

        ## If the echo level is high enough, print the complete list of settings used to run the simualtion
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def RunSolutionLoop(self):
        '''Break Solution Loop when Buckling Analysis is converged,
           or maximum step number is reached
        '''
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.step += 1
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
            if self._GetSolver().get_mechanical_solution_strategy().GetSolutionFoundFlag():
                break

    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop"""
        return self.step < self.nsteps

    def FinalizeSolutionStep(self):
        ''' This function is overriden to postprocess eigenvalues only every second load step.
        Print eigenvalues after every small load increment'''
        self._GetSolver().FinalizeSolutionStep()
        for process in self._GetListOfProcesses():
            if( process.__class__.__name__ != "PostprocessEigenvaluesProcess" ):
                process.ExecuteFinalizeSolutionStep()
            elif ( (self.step % 2 == 0) & (self.step > 0 ) ):
                process.ExecuteFinalizeSolutionStep()

    def OutputSolutionStep(self):
        ''' This function is overriden to print output only every second load step.
        Print output after every path following load step'''
        is_output_step = False
        for output_process in self._GetListOfOutputProcesses():
            if output_process.IsOutputStep():
                is_output_step = True
                break

        if is_output_step: # at least one of the output processes will print output
            for process in self._GetListOfProcesses():
                process.ExecuteBeforeOutputStep()

            for output_process in self._GetListOfOutputProcesses():
                if( output_process.IsOutputStep() & (self.step % 2 == 1) ):
                    output_process.PrintOutput()

            for process in self._GetListOfProcesses():
                process.ExecuteAfterOutputStep()

