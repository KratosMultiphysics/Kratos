from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

# Importing the solvers
import KratosMultiphysics.ExternalSolversApplication

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class PotentialFlowAnalysis(AnalysisStage):
    '''Main script for potential flow simulations.'''

    def __init__(self,model,project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        solver_settings = project_parameters["solver_settings"]

        if solver_settings.Has("domain_size") and project_parameters["problem_data"].Has("domain_size"):
            warn_msg  = '"domain_size" defined both in "problem_data" and "solver_settings"!'
            warn_msg += 'the definition in the "solver_settings" will be employed'
            KratosMultiphysics.Logger.PrintWarning("PotentialFlowAnalysis", warn_msg)

        if solver_settings.Has("model_part_name") and project_parameters["problem_data"].Has("model_part_name"):
            warn_msg  = '"model_part_name" defined both in "problem_data" and "solver_settings"!'
            warn_msg += 'the definition in the "solver_settings" will be employed'
            KratosMultiphysics.Logger.PrintWarning("PotentialFlowAnalysis", warn_msg)

        if not solver_settings.Has("domain_size") and project_parameters["problem_data"].Has("domain_size"):
            KratosMultiphysics.Logger.PrintWarning("PotentialFlowAnalysis", "Using the old way to pass the domain_size, this will be removed!")
            solver_settings.AddEmptyValue("domain_size")
            solver_settings["domain_size"].SetInt(project_parameters["problem_data"]["domain_size"].GetInt())

        if not solver_settings.Has("model_part_name") and project_parameters["problem_data"].Has("model_part_name"):
            KratosMultiphysics.Logger.PrintWarning("PotentialFlowAnalysis", "Using the old way to pass the model_part_name, this will be removed!")
            solver_settings.AddEmptyValue("model_part_name")
            solver_settings["model_part_name"].SetString(project_parameters["problem_data"]["model_part_name"].GetString())

        super(PotentialFlowAnalysis,self).__init__(model,project_parameters)

    def _CreateSolver(self):
        import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_solver as potential_flow_solver
        return potential_flow_solver.CreateSolver(self.model, self.project_parameters["solver_settings"])

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(PotentialFlowAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["initial_conditions_process_list", "boundary_conditions_process_list", "auxiliar_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                info_msg  = "Using the old way to create the processes, this will be removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                KratosMultiphysics.Logger.PrintWarning("PotentialFlowAnalysis", info_msg)
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
                info_msg  = "Using the old way to create the gid-output, this will be removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                KratosMultiphysics.Logger.PrintWarning("PotentialFlowAnalysis", info_msg)
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
                               self.project_parameters["problem_data"]["problem_name"].GetString(
        ),
            self.project_parameters["output_configuration"])

        return output

    def RunSolutionLoop(self):
        self.InitializeSolutionStep()
        self._GetSolver().Predict()
        self._GetSolver().SolveSolutionStep()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

