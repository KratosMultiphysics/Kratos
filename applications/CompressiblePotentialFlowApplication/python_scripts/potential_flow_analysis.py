# Importing Kratos
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class PotentialFlowAnalysis(AnalysisStage):
    '''Main script for potential flow simulations.'''

    def _CreateSolver(self):
        if self.project_parameters["solver_settings"]["solver_type"].GetString()=="potential_flow":
            import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_solver as potential_flow_solver
            return potential_flow_solver.CreateSolver(self.model, self.project_parameters["solver_settings"])
        elif self.project_parameters["solver_settings"]["solver_type"].GetString()=="ale_potential_flow":
            import KratosMultiphysics.CompressiblePotentialFlowApplication.ale_potential_flow_solver as ale_potential_flow_solver
            return ale_potential_flow_solver.CreateSolver(self.model, self.project_parameters["solver_settings"], self.parallel_type)
        elif self.project_parameters["solver_settings"]["solver_type"].GetString()=="adjoint_potential_flow":
            import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_adjoint_solver as adjoint_solver
            return adjoint_solver.CreateSolver(self.model, self.project_parameters["solver_settings"])
        else:
            raise Exception("Solver type '"+str(self.project_parameters["solver_settings"]["solver_type"].GetString())+"' not added. Please specify an available solver")


    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(PotentialFlowAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["boundary_conditions_process_list", "auxiliar_process_list"]
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
            from KratosMultiphysics.mpi.distributed_gid_output_process import DistributedGiDOutputProcess as OutputProcess

        output = OutputProcess(self._GetSolver().GetComputingModelPart(),
                               self.project_parameters["problem_data"]["problem_name"].GetString(
        ),
            self.project_parameters["output_configuration"])

        return output

    # def RunSolutionLoop(self):
    #     self.InitializeSolutionStep()
    #     self._GetSolver().Predict()
    #     self._GetSolver().SolveSolutionStep()
    #     self.FinalizeSolutionStep()
    #     self.OutputSolutionStep()

