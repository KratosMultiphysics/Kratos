# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM
import math

from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural as structural_solvers

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class AdjointStructuralMechanicsDynamicAnalysis(AnalysisStage):
    """
    Main script for adjoint sensitivity calculation in StructuralMechanicsApplication simulations.

    """

    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initializing the parameters
        solver_settings = project_parameters["solver_settings"]

        if solver_settings.Has("domain_size") and project_parameters["problem_data"].Has("domain_size"):
            raise Exception("AdjointStructuralMechanicsDynamicAnalysis: " + '"domain_size" defined both in "problem_data" and "solver_settings"!')

        if solver_settings.Has("model_part_name") and project_parameters["problem_data"].Has("model_part_name"):
            raise Exception("AdjointStructuralMechanicsDynamicAnalysis: " + '"model_part_name" defined both in problem_data" and "solver_settings"!')

        if solver_settings.Has("time_stepping") and project_parameters["problem_data"].Has("time_Step"):
            raise Exception("AdjointStructuralMechanicsDynamicAnalysis: " + '"time_stepping" defined both in "problem_data" and "solver_settings"!')

        if not solver_settings.Has("time_stepping"):
            raise Exception("AdjointStructuralMechanicsDynamicAnalysis: Using the old way to pass the time_step, this was removed!")

        if not solver_settings.Has("domain_size"):
            raise Exception("AdjointStructuralMechanicsDynamicAnalysis: Using the old way to pass the domain_size, this was removed!")

        # Detect is a contact problem
        # NOTE: We have a special treatment for contact problems due to the way the convergence info is printed (in a table). Not doing this will provoque that the table is discontinous (and not fancy and eye-candy)
        solver_settings = project_parameters["solver_settings"]
        self.contact_problem = solver_settings.Has("contact_settings") or solver_settings.Has("mpc_contact_settings")

        super().__init__(model, project_parameters)

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
        print("delta time is ", delta_time)
        problem_data = project_parameters["problem_data"]
        self.start_time = problem_data["start_time"].GetDouble() + delta_time
        self.end_time = problem_data["end_time"].GetDouble() + delta_time

        # now set start_time, end_time of the problem_data bt flipping them
        self.project_parameters["problem_data"]["start_time"].SetDouble(self.end_time)
        self.project_parameters["problem_data"]["end_time"].SetDouble(self.start_time) 

        self.num_of_measurements = math.ceil((self.end_time - self.start_time) / delta_time)
        print("num of measurements are ", self.num_of_measurements)
        if self.num_of_measurements <= 0:
            raise RuntimeError(f"Number of measurements are {self.num_of_measurements}. Make sure they are > 0.")

    def Initialize(self):
        """ Initializing the Analysis """
        super().Initialize()

        # In case of contact problem
        if self.contact_problem:
            self._GetSolver().SetEchoLevel(self.echo_level)
            # To avoid many prints
            if self.echo_level == 0:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

         # dummy time step to correctly calculate DELTA_TIME
        #self._GetSolver().main_model_part.CloneTimeStep(self.time)

    # def Run(self):
    #     self._GetSolver().main_model_part.ProcessInfo[Kratos.STEP] = self.num_of_measurements + 1
    #     super().Run()

    def RunSolutionLoop(self):
        self._GetSolver().main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.num_of_measurements + 1
        super().RunSolutionLoop()

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """

        # In case of contact problem
        if self.contact_problem:
            # First we check if one of the output processes will print output in this step this is done to save computation in case none of them will print
            is_output_step = False
            for output_process in self._GetListOfOutputProcesses():
                if output_process.IsOutputStep():
                    is_output_step = True
                    break

            if is_output_step:
                # Informing the output will be created
                KratosMultiphysics.Logger.PrintWarning("AdjointStructuralMechanicsDynamicAnalysis", "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
                KratosMultiphysics.Logger.PrintWarning("AdjointStructuralMechanicsDynamicAnalysis", "TIME: ", self.time)

        # Creating output
        super().OutputSolutionStep()


    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        Note that as the adjoint problem is solved backward in time, the stopping criteria is current time being larger than the start one
        """
        #print("Inside Advancining loop check ", self.time, " vs ", self.start_time)
        return self.time > self.start_time

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not already in the model) """
        ## Solver construction
        return structural_solvers.CreateSolver(self.model, self.project_parameters)

    # def _CreateProcesses(self, parameter_name, initialization_order):
    #     """Create a list of Processes
    #     This method is TEMPORARY to not break existing code
    #     It will be removed in the future
    #     """
    #     list_of_processes = super()._CreateProcesses(parameter_name, initialization_order)

    #     if parameter_name == "processes":
    #         processes_block_names = ["constraints_process_list", "loads_process_list", "list_other_processes", "json_output_process",
    #             "json_check_process", "check_analytic_results_process", "contact_process_list"]
    #         if len(list_of_processes) == 0: # Processes are given in the old format (or no processes are specified)
    #             for process_name in processes_block_names:
    #                 if self.project_parameters.Has(process_name):
    #                     info_msg  = "Using the old way to create the processes, this was removed!\n"
    #                     info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
    #                     info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
    #                     info_msg += "for a description of the new format"
    #                     raise Exception("AdjointStructuralMechanicsDynamicAnalysis: " + info_msg)

    #         else: # Processes are given in the new format
    #             for process_name in processes_block_names:
    #                 if self.project_parameters.Has(process_name):
    #                     raise Exception("Mixing of process initialization is not allowed!")
    #     elif parameter_name == "output_processes":
    #         if self.project_parameters.Has("output_configuration"):
    #             info_msg  = "Using the old way to create the gid-output, this was removed!\n"
    #             info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
    #             info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
    #             info_msg += "for a description of the new format"
    #             raise Exception("AdjointStructuralMechanicsDynamicAnalysis: " + info_msg)
    #     else:
    #         raise NameError("wrong parameter name")

    #     return list_of_processes

    def _GetSimulationName(self):
        return "::[KSM Adjoint Dynamic Simulation]:: "

if __name__ == "__main__":
    from sys import argv

    primal_parameter_file_name = None
    adjoint_parameter_file_name = None

    parameters = KratosMultiphysics.Parameters(r'''{}''')

    if len(argv) == 2:
        adjoint_parameter_file_name = argv[1]
    elif len(argv) == 3:
        primal_parameter_file_name = argv[1]
        adjoint_parameter_file_name = argv[2]
    else:
        err_msg =  'Unexpected amount of input arguments!\n'
        err_msg += 'To run the primal structural problem followed by the adjoint solution, provide both parameter files:\n'
        err_msg += '    "python3 adjoint_structural_mechanics_dynamic_analysis.py <primal-parameter-file>.json <adjoint-parameter-file>.json"\n'
        err_msg += 'To run only the adjoint problem, provide only the adjoint parameter file:\n'
        err_msg += '    "python3 adjoint_structural_mechanics_dynamic_analysis.py <adjoint-parameter-file>.json"\n'
        raise Exception(err_msg)

    if primal_parameter_file_name is not None:
        with open(primal_parameter_file_name,'r') as primal_parameter_file:
            parameters.AddValue("primal_settings", KratosMultiphysics.Parameters(primal_parameter_file.read()))
    else:
        parameters.AddEmptyValue("primal_settings")

    with open(adjoint_parameter_file_name,'r') as adjoint_parameter_file:
        parameters.AddValue("adjoint_settings", KratosMultiphysics.Parameters(adjoint_parameter_file.read()))

    model = KratosMultiphysics.Model()

    if primal_parameter_file_name is not None:
        primal_simulation = StructuralMechanicsAnalysis(model,parameters["primal_settings"])
        primal_simulation.Run()

    adjoint_model = KratosMultiphysics.Model()
    adjoint_simulation = AdjointStructuralMechanicsDynamicAnalysis(adjoint_model,parameters["adjoint_settings"])
    adjoint_simulation.Run()

