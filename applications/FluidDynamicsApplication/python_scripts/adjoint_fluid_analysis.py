from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
import KratosMultiphysics.AdjointFluidApplication as KAdjoint

from analysis_stage import AnalysisStage

class AdjointFluidAnalysis(AnalysisStage):
    '''Main script for adjoint sensitivity optimization in fluid dynamics simulations.'''

    def __init__(self,model,parameters):
        # Deprecation warnings
        solver_settings = parameters["solver_settings"]
        if not solver_settings.Has("domain_size"):
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Using the old way to pass the domain_size, this will be removed!")
            solver_settings.AddEmptyValue("domain_size")
            solver_settings["domain_size"].SetInt(parameters["problem_data"]["domain_size"].GetInt())

        if not solver_settings.Has("model_part_name"):
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Using the old way to pass the model_part_name, this will be removed!")
            solver_settings.AddEmptyValue("model_part_name")
            solver_settings["model_part_name"].SetString(parameters["problem_data"]["model_part_name"].GetString())
        
        if not parameters["problem_data"].Has("end_time"):
            parameters["problem_data"].AddEmptyValue("end_time")
            parameters["problem_data"]["end_time"].SetDouble( \
                            parameters["problem_data"]["start_step"].GetDouble() + \
                            parameters["problem_data"]["nsteps"].GetInt()*solver_settings["time_stepping"]["time_step"].GetDouble()
                        )

        if not parameters["problem_data"].Has("start_time"):
            parameters["problem_data"].AddEmptyValue("start_time")
            parameters["problem_data"]["start_time"].SetDouble( \
                            parameters["problem_data"]["start_step"].GetDouble() \
                            )
        self.number_of_steps = parameters["problem_data"]["nsteps"].GetInt()

        super(AdjointFluidAnalysis, self).__init__(model, parameters)

        # If this is an MPI run, load the distributed memory modules
        if (self.parallel_type == "MPI"):
            from KratosMultiphysics.mpi import mpi
            import KratosMultiphysics.MetisApplication
            import KratosMultiphysics.TrilinosApplication
            self.is_printing_rank = (mpi.rank == 0)
        else:
            self.is_printing_rank = True

    def _CreateSolver(self):
        import python_solvers_wrapper_adjoint_fluid
        return python_solvers_wrapper_adjoint_fluid.CreateSolver(self.model, self.project_parameters)            

    def RunSolutionLoop(self):
        """Note that the adjoint problem is solved in reverse time
        """

        # dummy time step to correctly calculate DELTA_TIME
        dt = self._solver._ComputeDeltaTime()
        new_time = self.time - dt
        self._solver.main_model_part.CloneTimeStep(new_time)

        for step in range(self.number_of_steps):
            self.time = self._solver.AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._solver.Predict()
            self._solver.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

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
