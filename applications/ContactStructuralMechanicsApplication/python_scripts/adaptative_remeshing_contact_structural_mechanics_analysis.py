from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Other imports
import sys

# Import the base structural analysis
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis as BaseClass

class AdaptativeRemeshingContactStructuralMechanicsAnalysis(BaseClass):
    """
    This class is the main-script of the ContactStructuralMechanicsApplication when using adaptative remeshing put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):

        # Construct the base analysis.
        default_params = KM.Parameters("""
        {
            "max_iteration" : 1,
            "analysis_type" : "linear",
            "contact_settings" :
            {
                "fancy_convergence_criterion" : false
            }
        }
        """)
        if project_parameters["solver_settings"].Has("max_iteration"):
            self.non_linear_iterations = project_parameters["solver_settings"]["max_iteration"].GetInt()
        else:
            self.non_linear_iterations = 10
            project_parameters["solver_settings"].AddValue("max_iteration", default_params["max_iteration"])
        if project_parameters["solver_settings"].Has("analysis_type"):
            project_parameters["solver_settings"]["analysis_type"].SetString("linear")
        else:
            project_parameters["solver_settings"].AddValue("analysis_type", default_params["analysis_type"])
        if project_parameters["solver_settings"].Has("contact_settings"):
            if project_parameters["solver_settings"]["contact_settings"].Has("fancy_convergence_criterion"):
                project_parameters["solver_settings"]["contact_settings"]["fancy_convergence_criterion"].SetBool(False)
            else:
                project_parameters["solver_settings"]["contact_settings"].AddValue("fancy_convergence_criterion", default_params["contact_settings"]["fancy_convergence_criterion"])
        else:
            project_parameters["solver_settings"].AddValue("contact_settings", default_params["contact_settings"])
        self.process_remesh = False
        if project_parameters.Has("recursive_remeshing_process"):
            self.process_remesh = True
        if project_parameters.Has("processes"):
            if project_parameters["processes"].Has("recursive_remeshing_process"):
                self.process_remesh = True
        if not self.process_remesh:
            project_parameters["solver_settings"]["analysis_type"].SetString("linear")
        super(AdaptativeRemeshingContactStructuralMechanicsAnalysis, self).__init__(model, project_parameters)

    def Initialize(self):
        """ Initializing the Analysis """
        super(AdaptativeRemeshingContactStructuralMechanicsAnalysis, self).Initialize()
        computing_model_part = self._GetSolver().GetComputingModelPart()
        if not self.process_remesh:
            convergence_criteria = self._GetSolver().get_convergence_criterion()
            convergence_criteria.Initialize(computing_model_part)

        # Ensuring to have conditions on the BC before remesh
        is_surface = False
        for elem in computing_model_part.Elements:
            geom = elem.GetGeometry()
            if geom.WorkingSpaceDimension() != geom.LocalSpaceDimension():
                is_surface = True
            break

        if not is_surface:
            list_model_parts = []
            # We need to detect the conditions in the boundary conditions
            if self.project_parameters.Has("constraints_process_list"):
                constraints_process_list = self.project_parameters["constraints_process_list"]
                for i in range(0,constraints_process_list.size()):
                    item = constraints_process_list[i]
                    list_model_parts.append(item["Parameters"]["model_part_name"].GetString())
            skin_detection_parameters = KM.Parameters("""
            {
                "list_model_parts_to_assign_conditions" : []
            }
            """)
            for name_mp in list_model_parts:
                skin_detection_parameters["list_model_parts_to_assign_conditions"].Append(name_mp)

            if computing_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2:
                detect_skin = KM.SkinDetectionProcess2D(computing_model_part, skin_detection_parameters)
            else:
                detect_skin = KM.SkinDetectionProcess3D(computing_model_part, skin_detection_parameters)
            detect_skin.Execute()
        self._GetSolver().SetEchoLevel(self.echo_level)

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        # If we remesh using a process
        computing_model_part = self._GetSolver().GetComputingModelPart()
        root_model_part = computing_model_part.GetRootModelPart()
        if self.process_remesh:
            while self.time < self.end_time:
                self.time = self._GetSolver().AdvanceInTime(self.time)
                # We reinitialize if remeshed previously
                if root_model_part.Is(KM.MODIFIED):
                    self._ReInitializeSolver()
                self.InitializeSolutionStep()
                # We reinitialize if remeshed on the InitializeSolutionStep
                if root_model_part.Is(KM.MODIFIED):
                    self._ReInitializeSolver()
                    self.InitializeSolutionStep()
                self._GetSolver().Predict()
                self._GetSolver().SolveSolutionStep()
                self._transfer_slave_to_master()
                self.FinalizeSolutionStep()
                self.OutputSolutionStep()
        else: # Remeshing adaptively
            metric_process = self._GetSolver().get_metric_process()
            remeshing_process = self._GetSolver().get_remeshing_process()
            convergence_criteria = self._GetSolver().get_convergence_criterion()
            builder_and_solver = self._GetSolver().get_builder_and_solver()
            mechanical_solution_strategy = self._GetSolver().get_mechanical_solution_strategy()

            while self.time < self.end_time:
                self.time = self._GetSolver().AdvanceInTime(self.time)
                non_linear_iteration = 1
                while non_linear_iteration <= self.non_linear_iterations:
                    if computing_model_part.Is(KM.MODIFIED):
                        self._GetSolver().Clear()
                        # WE RECOMPUTE THE PROCESSES AGAIN
                        # Processes initialization
                        for process in self._GetListOfProcesses():
                            process.ExecuteInitialize()
                        ## Processes before the loop
                        for process in self._GetListOfProcesses():
                            process.ExecuteBeforeSolutionLoop()
                        ## Processes of initialize the solution step
                        for process in self._GetListOfProcesses():
                            process.ExecuteInitializeSolutionStep()
                    if non_linear_iteration == 1 or computing_model_part.Is(KM.MODIFIED):
                        self.InitializeSolutionStep()
                        self._GetSolver().Predict()
                        computing_model_part.Set(KM.MODIFIED, False)
                    computing_model_part.ProcessInfo.SetValue(KM.NL_ITERATION_NUMBER, non_linear_iteration)
                    is_converged = convergence_criteria.PreCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                    self._GetSolver().SolveSolutionStep()
                    is_converged = convergence_criteria.PostCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                    self._transfer_slave_to_master()
                    self.FinalizeSolutionStep()
                    if is_converged:
                        KM.Logger.PrintInfo(self._GetSimulationName(), "Adaptative strategy converged in ", non_linear_iteration, "iterations" )
                        break
                    elif non_linear_iteration == self.non_linear_iterations:
                        KM.Logger.PrintInfo(self._GetSimulationName(), "Adaptative strategy not converged after ", non_linear_iteration, "iterations" )
                        break
                    else:
                        # Before remesh we set the flag INTERFACE to the conditions (we need edges to preserve submodelparts)
                        KM.VariableUtils().SetFlag(KM.INTERFACE, True, computing_model_part.GetSubModelPart("Contact").Conditions)
                        # We remove the contact model part to avoid problems (it will  be recomputed later)
                        contact_model_part = computing_model_part.GetSubModelPart("Contact")
                        for model_part in contact_model_part.SubModelParts:
                            contact_model_part.RemoveSubModelPart(model_part.Name)
                        computing_model_part.RemoveSubModelPart("ComputingContact")
                        # Now we can remesh
                        metric_process.Execute()
                        remeshing_process.Execute()
                        # We remove the contact model part to avoid problems (it will  be recomputed later)
                        computing_model_part.RemoveSubModelPart("Contact")
                        # Now we set as modified
                        computing_model_part.Set(KM.MODIFIED, True)
                        non_linear_iteration += 1
                self.OutputSolutionStep()

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

        # To avoid many prints
        if self.echo_level == 0:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

        ## Solver construction
        import python_solvers_wrapper_adaptative_remeshing_contact_structural
        return python_solvers_wrapper_adaptative_remeshing_contact_structural.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(AdaptativeRemeshingContactStructuralMechanicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["recursive_remeshing_process"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                KM.Logger.PrintWarning("AdaptativeRemeshingContactStructuralMechanicsAnalysis", "Using the old way to create the processes, this will be removed!")
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
            pass # Already added
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _transfer_slave_to_master(self):
        """ This method to transfer information from the slave side to the master side

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        computing_model_part = self._GetSolver().GetComputingModelPart()
        map_parameters = KM.Parameters("""
        {
            "echo_level"                       : 0,
            "destination_variable_historical"  : false,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "origin_variable"                  : "AUGMENTED_NORMAL_CONTACT_PRESSURE",
            "destination_variable"             : "CONTACT_PRESSURE",
            "integration_order"                : 2
        }
        """)

        interface_model_part = computing_model_part.GetSubModelPart("Contact")
        main_model_part = computing_model_part.GetRootModelPart()
        main_model_part.RemoveSubModelPart("ComputingContact")

        if interface_model_part.HasSubModelPart("SlaveSubModelPart"):
            slave_interface_model_part = interface_model_part.GetSubModelPart("SlaveSubModelPart")
        else:
            slave_interface_model_part = interface_model_part.CreateSubModelPart("SlaveSubModelPart")
            KM.FastTransferBetweenModelPartsProcess(slave_interface_model_part, interface_model_part, KM.FastTransferBetweenModelPartsProcess.EntityTransfered.NODES, KM.SLAVE).Execute()
            KM.FastTransferBetweenModelPartsProcess(slave_interface_model_part, interface_model_part, KM.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS, KM.SLAVE).Execute()
        if interface_model_part.HasSubModelPart("MasterSubModelPart"):
            master_interface_model_part = interface_model_part.GetSubModelPart("MasterSubModelPart")
        else:
            master_interface_model_part = interface_model_part.CreateSubModelPart("MasterSubModelPart")
            KM.FastTransferBetweenModelPartsProcess(master_interface_model_part, interface_model_part, KM.FastTransferBetweenModelPartsProcess.EntityTransfered.NODES, KM.MASTER).Execute()
            KM.FastTransferBetweenModelPartsProcess(master_interface_model_part, interface_model_part, KM.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS, KM.MASTER).Execute()

        mortar_mapping = KM.SimpleMortarMapperProcess(slave_interface_model_part, master_interface_model_part, map_parameters)
        mortar_mapping.Execute()

    def _ReInitializeSolver(self):
        """ This reinitializes after remesh """
        self._GetSolver().Clear()
        # WE INITIALIZE THE SOLVER
        self._GetSolver().Initialize()
        # WE RECOMPUTE THE PROCESSES AGAIN
        ## Processes initialization
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()
        ## Processes before the loop
        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()
        ## Processes of initialize the solution step
        for process in self._GetListOfProcesses():
            process.ExecuteInitializeSolutionStep()

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 adaptative_contact_structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 adaptative_contact_structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    AdaptativeRemeshingContactStructuralMechanicsAnalysis(project_parameters_file_name).Run()
