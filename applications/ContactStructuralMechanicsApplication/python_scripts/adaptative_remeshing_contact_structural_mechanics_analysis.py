from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM

# Other imports
import sys

# Import the base structural analysis
from KratosMultiphysics.ContactStructuralMechanicsApplication.contact_structural_mechanics_analysis import ContactStructuralMechanicsAnalysis as BaseClass

# Import auxiliar methods
from KratosMultiphysics.auxiliar_methods_adaptative_remeshing import AuxiliarMethodsAdaptiveRemeshing

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
        if project_parameters.Has("mesh_adaptivity_processes"):
            self.process_remesh = True
        if project_parameters.Has("processes"):
            if project_parameters["processes"].Has("mesh_adaptivity_processes"):
                self.process_remesh = True
        if not self.process_remesh:
            project_parameters["solver_settings"]["analysis_type"].SetString("linear")
        super(AdaptativeRemeshingContactStructuralMechanicsAnalysis, self).__init__(model, project_parameters)

        # Create utilities
        self.adaptive_utilities = AuxiliarMethodsAdaptiveRemeshing(self)

    def Initialize(self):
        """ Initializing the Analysis """
        super(AdaptativeRemeshingContactStructuralMechanicsAnalysis, self).Initialize()
        self.adaptive_utilities.AdaptativeRemeshingDetectBoundary()

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        # If we remesh using a process
        if self.process_remesh:
            self.adaptive_utilities.AdaptativeRemeshingRunSolutionLoop()
        else: # Remeshing adaptively
            self.adaptive_utilities.SPRAdaptativeRemeshingRunSolutionLoop(True)

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
            processes_block_names = ["mesh_adaptivity_processes"]
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

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 adaptative_remeshing_contact_structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 adaptative_remeshing_contact_structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    AdaptativeRemeshingContactStructuralMechanicsAnalysis(project_parameters_file_name).Run()
