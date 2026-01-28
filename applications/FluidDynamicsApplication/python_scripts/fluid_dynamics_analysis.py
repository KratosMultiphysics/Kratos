
from sys import argv

import KratosMultiphysics as Kratos

from KratosMultiphysics.analysis_stage_with_solver import AnalysisStageWithSolver
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid

class FluidDynamicsAnalysis(AnalysisStageWithSolver):
    '''Main script for fluid dynamics simulations using the navier_stokes family of python solvers.'''

    def __init__(self,model,parameters):
        # Deprecation warnings
        solver_settings = parameters["solver_settings"]

        if solver_settings.Has("domain_size") and parameters["problem_data"].Has("domain_size"):
            raise Exception('FluidDynamicsAnalysis: "domain_size" defined both in "problem_data" and "solver_settings"!')

        if solver_settings.Has("model_part_name") and parameters["problem_data"].Has("model_part_name"):
            raise Exception('FluidDynamicsAnalysis: "model_part_name" defined both in "problem_data" and "solver_settings"!')

        if not solver_settings.Has("domain_size") and parameters["problem_data"].Has("domain_size"):
            raise Exception("FluidDynamicsAnalysis: Using the old way to pass the domain_size, this was removed!")

        if not solver_settings.Has("model_part_name") and parameters["problem_data"].Has("model_part_name"):
            raise Exception("FluidDynamicsAnalysis: Using the old way to pass the model_part_name, this was removed!")

        super(FluidDynamicsAnalysis,self).__init__(model,parameters)

    def _CreateSolver(self):
        return python_solvers_wrapper_fluid.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(FluidDynamicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["gravity", "initial_conditions_process_list", "boundary_conditions_process_list", "auxiliar_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format (or no processes are specified)
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        info_msg  = "Using the old way to create the processes, this was removed!\n"
                        info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                        info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                        info_msg += "for a description of the new format"
                        raise Exception("FluidDynamicsAnalysis: " + info_msg)
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not allowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                info_msg  = "Using the old way to create the gid-output, this was removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                raise Exception("FluidDynamicsAnalysis: " + info_msg)
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _GetOrderOfProcessesInitialization(self):
        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        return ["gravity",
                "initial_conditions_process_list",
                "boundary_conditions_process_list",
                "auxiliar_process_list"]

    def _GetSimulationName(self):
        return "Fluid Dynamics Analysis"

if __name__ == '__main__':
    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = FluidDynamicsAnalysis(model,parameters)
    simulation.Run()
