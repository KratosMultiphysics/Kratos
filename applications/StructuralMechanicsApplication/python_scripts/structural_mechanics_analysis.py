from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Importing the solvers (if available)
try:
    import KratosMultiphysics.ExternalSolversApplication
    KratosMultiphysics.Logger.PrintInfo("ExternalSolversApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("ExternalSolversApplication", "not imported")
try:
    import KratosMultiphysics.EigenSolversApplication
    KratosMultiphysics.Logger.PrintInfo("EigenSolversApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("EigenSolversApplication", "not imported")

# Importing the base class
from analysis_stage import AnalysisStage

# Other imports
import sys

class StructuralMechanicsAnalysis(AnalysisStage):
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # In case of the old format, where dT is given under "problem_data", copy it to solver settings
        solver_settings = project_parameters["solver_settings"]
        if not solver_settings.Has("time_stepping"):
            KratosMultiphysics.Logger.PrintInfo("StructuralMechanicsAnalysis", "Using the old way to pass the time_step, this will be removed!")
            time_stepping_params = KratosMultiphysics.Parameters("{}")
            time_stepping_params.AddValue("time_step", project_parameters["problem_data"]["time_step"])
            solver_settings.AddValue("time_stepping", time_stepping_params)

        # Create the ModelPart
        # Note that this in temporary and will be done through the model in the future
        main_model_part_name = project_parameters["problem_data"]["model_part_name"].GetString()
        if model.HasModelPart(main_model_part_name):
            self.main_model_part = model[main_model_part_name]
        else:
            self.main_model_part = KratosMultiphysics.ModelPart(main_model_part_name)
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,
                                                      project_parameters["problem_data"]["domain_size"].GetInt())

        super(StructuralMechanicsAnalysis, self).__init__(model, project_parameters)


    def OutputSolutionStep(self):
        super(StructuralMechanicsAnalysis, self).OutputSolutionStep()

        self._GetSolver().SaveRestart()


    #### Internal functions ####
    def _CreateSolver(self, external_model_part=None):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        import python_solvers_wrapper_structural
        return python_solvers_wrapper_structural.CreateSolver(self.main_model_part, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is temporary to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(StructuralMechanicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        if parameter_name == "processes":
            old_processes_names = ["constraints_process_list", "loads_process_list", "list_other_processes", "json_output_process"
                "json_check_process", "check_analytic_results_process", "contact_process_list"]
            if len(list_of_processes) == 0:
                KratosMultiphysics.Logger.PrintInfo("StructuralMechanicsAnalysis", "Using the old way to create the processes, this will be removed!")
                from process_factory import KratosProcessFactory
                factory = KratosProcessFactory(self.model)
                for process_name in old_processes_names:
                    if (self.project_parameters.Has(process_name) is True):
                        list_of_processes += factory.ConstructListOfProcesses(self.project_parameters[process_name])
            else:
                for process_name in old_processes_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not alowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                #KratosMultiphysics.Logger.PrintInfo("StructuralMechanicsAnalysis", "Using the old way to create the gid-output, this will be removed!")
                gid_output= self._SetUpGiDOutput()
                list_of_processes += [gid_output,]
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _SetUpGiDOutput(self):
        '''Initialize a GiD output instance'''
        self.__CheckForDeprecatedGiDSettings()
        if self.parallel_type == "OpenMP":
            from gid_output_process import GiDOutputProcess as OutputProcess
        elif self.parallel_type == "MPI":
            from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

        gid_output = OutputProcess(self._GetSolver().GetComputingModelPart(),
                                   self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                   self.project_parameters["output_configuration"])

        return gid_output

    def _GetSimulationName(self):
        return "::[KSM Simulation]:: "

    def __CheckForDeprecatedGiDSettings(self):
        if self.project_parameters["output_configuration"].Has("result_file_configuration"):
            res_file_config = self.project_parameters["output_configuration"]["result_file_configuration"]
            if res_file_config.Has("nodal_results"):
                nodal_res = res_file_config["nodal_results"]
                for i in range(nodal_res.size()):
                    var_name = nodal_res[i].GetString()
                    if var_name == "TORQUE":
                        err_msg  = 'Requesting output for "TORQUE" which is not available any more\n'
                        err_msg += 'It was renamed to "MOMENT_REACTION"'
                        raise Exception(err_msg)

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Run()
