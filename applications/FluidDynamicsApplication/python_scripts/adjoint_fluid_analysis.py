from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KFluid
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass

from analysis_stage import AnalysisStage
from fluid_dynamics_analysis import FluidDynamicsAnalysis

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

        self.is_printing_rank = True
        ## Import parallel modules if needed
        if (parameters["problem_data"]["parallel_type"].GetString() == "MPI"):
            from KratosMultiphysics.mpi import mpi
            import KratosMultiphysics.MetisApplication as MetisApplication
            import KratosMultiphysics.TrilinosApplication as TrilinosApplication
            self.is_printing_rank = (mpi.rank == 0)

        super(AdjointFluidAnalysis, self).__init__(model, parameters)

    def Initialize(self):
        super(AdjointFluidAnalysis, self).Initialize()

        # dummy time step to correctly calculate DELTA_TIME
        self._GetSolver().main_model_part.CloneTimeStep(self.time)

    def _CreateSolver(self):
        import python_solvers_wrapper_adjoint_fluid
        return python_solvers_wrapper_adjoint_fluid.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(AdjointFluidAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        if parameter_name == "processes":
            processes_block_names = ["gravity", "initial_conditions_process_list", "boundary_conditions_process_list", "auxiliar_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Using the old way to create the processes, this will be removed!")
                from process_factory import KratosProcessFactory
                factory = KratosProcessFactory(self.model)
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        list_of_processes += factory.ConstructListOfProcesses(self.project_parameters[process_name])
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not alowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                #KratosMultiphysics.Logger.PrintInfo("FluidDynamicsAnalysis", "Using the old way to create the gid-output, this will be removed!")
                gid_output= self._SetUpGiDOutput()
                list_of_processes += [gid_output,]
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _SetUpGiDOutput(self):
        '''Initialize a GiD output instance'''
        if self.parallel_type == "OpenMP":
            from gid_output_process import GiDOutputProcess as OutputProcess
        elif self.parallel_type == "MPI":
            from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

        output = OutputProcess(self._GetSolver().GetComputingModelPart(),
                                self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                self.project_parameters["output_configuration"])

        return output

    def RunSolutionLoop(self):
        """Note that the adjoint problem is solved in reverse time
        """
        for step in range(self.number_of_steps):
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
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
