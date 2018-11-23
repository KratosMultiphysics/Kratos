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
try:
    import KratosMultiphysics.MeshingApplication
    KratosMultiphysics.Logger.PrintInfo("MeshingApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("MeshingApplication", "not imported")

# Importing the base class
from analysis_stage import AnalysisStage

class StructuralMechanicsAnalysis(AnalysisStage):
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        solver_settings = project_parameters["solver_settings"]

        if solver_settings.Has("domain_size") and project_parameters["problem_data"].Has("domain_size"):
            warn_msg  = '"domain_size" defined both in "problem_data" and "solver_settings"!'
            warn_msg += 'the definition in the "solver_settings" will be employed'
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", warn_msg)

        if solver_settings.Has("model_part_name") and project_parameters["problem_data"].Has("model_part_name"):
            warn_msg  = '"model_part_name" defined both in problem_data" and "solver_settings"!'
            warn_msg += 'the definition in the "solver_sett"ings" will be employed'
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", warn_msg)

        if solver_settings.Has("time_stepping") and project_parameters["problem_data"].Has("time_Step"):
            warn_msg  = '"time_stepping" defined both in "problem_data" and "solver_settings"!'
            warn_msg += 'the definition in the "solver_settings" will be employed'
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", warn_msg)

        if not solver_settings.Has("time_stepping"):
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", "Using the old way to pass the time_step, this will be removed!")
            time_stepping_params = KratosMultiphysics.Parameters("{}")
            time_stepping_params.AddValue("time_step", project_parameters["problem_data"]["time_step"])
            solver_settings.AddValue("time_stepping", time_stepping_params)

        if not solver_settings.Has("domain_size"):
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", "Using the old way to pass the domain_size, this will be removed!")
            solver_settings.AddEmptyValue("domain_size")
            solver_settings["domain_size"].SetInt(project_parameters["problem_data"]["domain_size"].GetInt())

        if not solver_settings.Has("model_part_name"):
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", "Using the old way to pass the model_part_name, this will be removed!")
            solver_settings.AddEmptyValue("model_part_name")
            solver_settings["model_part_name"].SetString(project_parameters["problem_data"]["model_part_name"].GetString())

        # Import parallel modules if needed
        # has to be done before the base-class constuctor is called (in which the solver is constructed)
        if (project_parameters["problem_data"]["parallel_type"].GetString() == "MPI"):
            import KratosMultiphysics.MetisApplication as MetisApplication
            import KratosMultiphysics.TrilinosApplication as TrilinosApplication

        super(StructuralMechanicsAnalysis, self).__init__(model, project_parameters)


    def Check(self):
        super(StructuralMechanicsAnalysis, self).Check()

        # performing some checks if the submodelparts used for the processes and
        # the material-assignments are being added to the ComputingModelPart
        solver_settings = self.project_parameters["solver_settings"]
        main_model_part_name = solver_settings["model_part_name"].GetString()

        # Checking if the material-submodelparts are added to the ComputingModelPart
        materials_filename = solver_settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""): # Materials are specified through a file
            # creating a list with the names of smps that will be added to the ComputingModelPart
            # note that the names here are WITHOUT the MainModelPart-Name
            domain_smp_param = solver_settings["problem_domain_sub_model_part_list"]
            list_domain_smp_names = [domain_smp_param[i].GetString() for i in range(domain_smp_param.size())]

            with open(materials_filename,'r') as materials_file: # reading the materials-file
                materials = KratosMultiphysics.Parameters(materials_file.read())

            for i in range(materials["properties"].size()):
                model_part_name = materials["properties"][i]["model_part_name"].GetString()
                if model_part_name.startswith(main_model_part_name): # removing the MainModelPart-Name
                    model_part_name = model_part_name.replace(main_model_part_name+".", "")
                if model_part_name not in list_domain_smp_names:
                    warn_msg  = 'The SubModelPart with name "' + model_part_name + '"\n'
                    warn_msg += 'is used for assigning materials but is not added to the ComputingModelPart!\n'
                    warn_msg += 'This can be done by adding it to "problem_domain_sub_model_part_list" '
                    warn_msg += 'in "solver_settings"\n'
                    KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis; Warning", warn_msg)

        if not self.project_parameters.Has("processes"): # TODO this check is only for backwards-compatibility => to be removed
            return

        # Checking if the processes-submodelparts are added to the ComputingModelPart
        # creating a list with the names of smps that will be added to the ComputingModelPart
        # note that the names here are WITHOUT the MainModelPart-Name
        processes_smp_param = solver_settings["processes_sub_model_part_list"]
        list_proc_smp_names = [processes_smp_param[i].GetString() for i in range(processes_smp_param.size())]

        for processes_block in self.project_parameters["processes"].values():
            for i_proc in range(processes_block.size()):
                process_params = processes_block[i_proc]["Parameters"]
                if process_params.Has("model_part_name"):
                    model_part_name = process_params["model_part_name"].GetString()
                    if model_part_name == main_model_part_name:
                        continue # skip the check for the MainModelPart
                    if model_part_name.startswith(main_model_part_name): # removing the MainModelPart-Name
                        model_part_name = model_part_name.replace(main_model_part_name+".", "")
                        KratosMultiphysics.Logger.PrintWarning(model_part_name)
                    if model_part_name not in list_proc_smp_names:
                        proc_name = processes_block[i_proc]["python_module"].GetString()
                        warn_msg  = 'The SubModelPart with name "' + model_part_name + '"\n'
                        warn_msg += 'is used for a process ("{}") \nbut is not added to the '.format(proc_name)
                        warn_msg += 'ComputingModelPart!\nThis can be done by adding it to '
                        warn_msg += '"processes_sub_model_part_list" in "solver_settings"\n'
                        KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis; Warning", warn_msg)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        import python_solvers_wrapper_structural
        return python_solvers_wrapper_structural.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(StructuralMechanicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["constraints_process_list", "loads_process_list", "list_other_processes", "json_output_process",
                "json_check_process", "check_analytic_results_process", "contact_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                info_msg  = "Using the old way to create the processes, this will be removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", info_msg)
                from process_factory import KratosProcessFactory
                factory = KratosProcessFactory(self.model)
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        list_of_processes += factory.ConstructListOfProcesses(self.project_parameters[process_name])
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not allowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                info_msg  = "Using the old way to create the gid-output, this will be removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                KratosMultiphysics.Logger.PrintInfo("StructuralMechanicsAnalysis", info_msg)
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
                        err_msg += 'It was renamed to "REACTION_MOMENT"'
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

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Run()
