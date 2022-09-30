import KratosMultiphysics as Kratos

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis
from KratosMultiphysics.RANSApplication.adjoint_rans_solver import AdjointRANSSolver

class AdjointRANSAnalysis(AnalysisStage):
    '''Main script for adjoint sensitivity optimization in fluid dynamics simulations.'''

    def __init__(self,model,parameters):
        super().__init__(model, parameters)

        problem_data = self.project_parameters["problem_data"]

        # compute delta time
        delta_time = self._GetSolver()._ComputeDeltaTime()
        if (delta_time >= 0.0):
            raise Exception("Adjoints are solved in reverse in time. Hence, delta_time should be negative. [ delta_time = " + str(delta_time) + "].")
        delta_time *= -1.0

        # update start time and end time. Adjoints are solved reverse in time,
        # hence we include an additional time step at the end
        # to properly initialize initial conditions for adjoint problem
        # and start time is also shifted by one to keep the same number
        # of steps
        self.start_time = problem_data["start_time"].GetDouble() + delta_time
        self.end_time = problem_data["end_time"].GetDouble() + delta_time

        # in the case of steady problems, is is only required to do calculations
        # in the last time step. therefore,
        if (self._GetSolver().is_steady):
            self.start_time = self.end_time - delta_time

        # now set start_time, end_time of the problem_data bt flipping them
        self.project_parameters["problem_data"]["start_time"].SetDouble(self.end_time)
        self.project_parameters["problem_data"]["end_time"].SetDouble(self.start_time)

    def Initialize(self):
        super().Initialize()

        # dummy time step to correctly calculate DELTA_TIME
        # and to initialize zero values as initial conditions for
        # adjoint variables
        self._GetSolver().main_model_part.CloneTimeStep(self.time)

    def _CreateSolver(self):
        return AdjointRANSSolver(self.model, self.project_parameters["solver_settings"])

    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        return self.time > self.start_time

    def _GetOrderOfProcessesInitialization(self):
        return ["initial_conditions_process_list",
                "boundary_conditions_process_list",
                "auxiliar_process_list"]

    def _GetSimulationName(self):
        return self.__class__.__name__

if __name__ == '__main__':
    from sys import argv

    primal_parameter_file_name = None
    adjoint_parameter_file_name = None

    parameters = Kratos.Parameters(r'''{}''')

    if len(argv) == 2:
        adjoint_parameter_file_name = argv[1]
    elif len(argv) == 3:
        primal_parameter_file_name = argv[1]
        adjoint_parameter_file_name = argv[2]
    else:
        err_msg =  'Unexpected amount of input arguments!\n'
        err_msg += 'To run the primal fluid problem followed by the adjoint solution, provide both parameter files:\n'
        err_msg += '    "python adjoint_fluid_analysis.py <primal-parameter-file>.json <adjoint-parameter-file>.json"\n'
        err_msg += 'To run only the adjoint problem, provide only the adjoint parameter file:\n'
        err_msg += '    "python adjoint_fluid_analysis.py <adjoint-parameter-file>.json"\n'
        raise Exception(err_msg)

    if primal_parameter_file_name is not None:
        with open(primal_parameter_file_name,'r') as primal_parameter_file:
            parameters.AddValue("primal_settings",Kratos.Parameters(primal_parameter_file.read()))
    else:
        parameters.AddEmptyValue("primal_settings")

    with open(adjoint_parameter_file_name,'r') as adjoint_parameter_file:
        parameters.AddValue("adjoint_settings", Kratos.Parameters(adjoint_parameter_file.read()))

    model = Kratos.Model()

    if primal_parameter_file_name is not None:
        primal_simulation = RANSAnalysis(model, parameters["primal_settings"])
        primal_simulation.Run()

    adjoint_model = Kratos.Model()
    adjoint_simulation = AdjointRANSAnalysis(adjoint_model, parameters["adjoint_settings"])
    adjoint_simulation.Run()
