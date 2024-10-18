
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KFluid

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_adjoint_fluid
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class AdjointFluidAnalysis(AnalysisStage):
    '''Main script for adjoint sensitivity optimization in fluid dynamics simulations.'''

    def __init__(self,model,parameters):
        super().__init__(model, parameters)

        problem_data = self.project_parameters["problem_data"]

        # check that old input is not provided
        if problem_data.Has("nsteps"):
            raise Exception("'nsteps' is no longer supported. Use 'start_time' and 'end_time' instead")
        if problem_data.Has("start_step"):
            raise Exception("'start_step' is no longer supported. Use 'start_time' and 'end_time' instead")
        if problem_data.Has("end_step"):
            raise Exception("'end_step' is no longer supported. Use 'start_time' and 'end_time' instead")

        # compute delta time
        delta_time = self._GetSolver()._ComputeDeltaTime()
        if delta_time >= 0.0:
            raise Exception("Adjoints are solved in reverse in time. Hence, delta_time should be negative. [ delta_time = " + str(delta_time) + "].")

        # update start time and end time. Adjoints are solved reverse in time,
        # hence we include an additional time step at the end
        # to properly initialize initial conditions for adjoint problem
        # and start time is also shifted by one to keep the same number
        # of steps
        delta_time *= -1.0
        self.start_time = problem_data["start_time"].GetDouble() + delta_time
        self.end_time = problem_data["end_time"].GetDouble() + delta_time

        # now set start_time, end_time of the problem_data bt flipping them
        self.project_parameters["problem_data"]["start_time"].SetDouble(self.end_time)
        self.project_parameters["problem_data"]["end_time"].SetDouble(self.start_time)

    def Initialize(self):
        super(AdjointFluidAnalysis, self).Initialize()

        # dummy time step to correctly calculate DELTA_TIME
        self._GetSolver().main_model_part.CloneTimeStep(self.time)

    def _CreateSolver(self):
        return python_solvers_wrapper_adjoint_fluid.CreateSolver(self.model, self.project_parameters)

    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        Note that as the adjoint problem is solved backward in time, the stopping criteria is current time being larger than the start one
        """
        return self.time > self.start_time

    def _GetOrderOfProcessesInitialization(self):
        return ["gravity",
                "initial_conditions_process_list",
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
        primal_simulation = FluidDynamicsAnalysis(model,parameters["primal_settings"])
        primal_simulation.Run()

    adjoint_model = Kratos.Model()
    adjoint_simulation = AdjointFluidAnalysis(adjoint_model,parameters["adjoint_settings"])
    adjoint_simulation.Run()
