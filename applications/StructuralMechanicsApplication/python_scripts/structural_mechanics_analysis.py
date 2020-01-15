from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural as structural_solvers

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class StructuralMechanicsAnalysis(AnalysisStage):
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        solver_settings = project_parameters["solver_settings"]

        if solver_settings.Has("domain_size") and project_parameters["problem_data"].Has("domain_size"):
            raise Exception("StructuralMechanicsAnalysis: " + '"domain_size" defined both in "problem_data" and "solver_settings"!')

        if solver_settings.Has("model_part_name") and project_parameters["problem_data"].Has("model_part_name"):
            raise Exception("StructuralMechanicsAnalysis: " + '"model_part_name" defined both in problem_data" and "solver_settings"!')

        if solver_settings.Has("time_stepping") and project_parameters["problem_data"].Has("time_Step"):
            raise Exception("StructuralMechanicsAnalysis: " + '"time_stepping" defined both in "problem_data" and "solver_settings"!')

        if not solver_settings.Has("time_stepping"):
            raise Exception("StructuralMechanicsAnalysis: Using the old way to pass the time_step, this was removed!")

        if not solver_settings.Has("domain_size"):
            raise Exception("StructuralMechanicsAnalysis: Using the old way to pass the domain_size, this was removed!")

        # Detect is a contact problem
        # NOTE: We have a special treatment for contact problems due to the way the convergence info is printed (in a table). Not doing this will provoque that the table is discontinous (and not fancy and eye-candy)
        solver_settings = project_parameters["solver_settings"]
        self.contact_problem = solver_settings.Has("contact_settings") or solver_settings.Has("mpc_contact_settings")

        if self.contact_problem:
            if solver_settings.Has("use_computing_model_part"):
                if not solver_settings["use_computing_model_part"].GetBool():
                    KM.Logger.PrintInfo("StructuralMechanicsAnalysis", 'For a contact problem the "ComputingModelPart" has to be used for now! Switching to True')
                    solver_settings["use_computing_model_part"].SetBool(True)
            else:
                solver_settings.AddEmptyValue("use_computing_model_part").SetBool(True)


        super(StructuralMechanicsAnalysis, self).__init__(model, project_parameters)

    def Initialize(self):
        """ Initializing the Analysis """
        super(StructuralMechanicsAnalysis, self).Initialize()

        # In case of contact problem
        if self.contact_problem:
            self._GetSolver().SetEchoLevel(self.echo_level)
            # To avoid many prints
            if self.echo_level == 0:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

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
                KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
                KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", "TIME: ", self.time)

        # Creating output
        super(StructuralMechanicsAnalysis, self).OutputSolutionStep()


    def Check(self):
        super(StructuralMechanicsAnalysis, self).Check()

        # performing some checks if the submodelparts used for the processes and
        # the material-assignments are being added to the ComputingModelPart
        solver_settings = self.project_parameters["solver_settings"]
        if not solver_settings["use_computing_model_part"].GetBool():
            return # no computing model part used, hence checks are not necessary

        main_model_part_name = solver_settings["model_part_name"].GetString()

        # Checking if the material-submodelparts are added to the ComputingModelPart
        materials_filename = solver_settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""): # Materials are specified through a file
            # creating a list with the names of smps that will be added to the ComputingModelPart
            # note that the names here are WITHOUT the MainModelPart-Name
            domain_smp_param = solver_settings["problem_domain_sub_model_part_list"]
            list_domain_mp_names = [domain_smp_param[i].GetString() for i in range(domain_smp_param.size())]

            if not main_model_part_name in list_domain_mp_names:
                # if the mainmodelpart is added to the computingmodelpart, then also all
                # submodelparts are added, no need to further check the submodelparts

                with open(materials_filename,'r') as materials_file: # reading the materials-file
                    materials = KratosMultiphysics.Parameters(materials_file.read())

                for i in range(materials["properties"].size()):
                    model_part_name = materials["properties"][i]["model_part_name"].GetString()
                    if model_part_name.startswith(main_model_part_name): # removing the MainModelPart-Name
                        model_part_name = model_part_name.replace(main_model_part_name+".", "")
                    if model_part_name not in list_domain_mp_names:
                        warn_msg  = 'The ModelPart with name "' + model_part_name + '"\n'
                        warn_msg += 'is used for assigning materials but is not added to the ComputingModelPart!\n'
                        warn_msg += 'This can be done by adding it to "problem_domain_sub_model_part_list" '
                        warn_msg += 'in "solver_settings"\n'
                        KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis; Warning", warn_msg)

        # Checking if the processes-submodelparts are added to the ComputingModelPart
        # creating a list with the names of smps that will be added to the ComputingModelPart
        # note that the names here are WITHOUT the MainModelPart-Name
        processes_smp_param = solver_settings["processes_sub_model_part_list"]
        list_proc_mp_names = [processes_smp_param[i].GetString() for i in range(processes_smp_param.size())]

        if not main_model_part_name in list_proc_mp_names:
            # if the mainmodelpart is added to the computingmodelpart, then also all
            # submodelparts are added, no need to further check the submodelparts

            for processes_block in self.project_parameters["processes"].values():
                for i_proc in range(processes_block.size()):
                    process_params = processes_block[i_proc]["Parameters"]
                    if process_params.Has("model_part_name"):
                        model_part_name = process_params["model_part_name"].GetString()
                        if model_part_name.startswith(main_model_part_name): # removing the MainModelPart-Name
                            model_part_name = model_part_name.replace(main_model_part_name+".", "")
                        if model_part_name not in list_proc_mp_names:
                            proc_name = processes_block[i_proc]["python_module"].GetString()
                            warn_msg  = 'The ModelPart with name "' + model_part_name + '"\n'
                            warn_msg += 'is used for a process ("{}") \nbut is not added to the '.format(proc_name)
                            warn_msg += 'ComputingModelPart!\nThis can be done by adding it to '
                            warn_msg += '"processes_sub_model_part_list" in "solver_settings"\n'
                            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis; Warning", warn_msg)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return structural_solvers.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(StructuralMechanicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["constraints_process_list", "loads_process_list", "list_other_processes", "json_output_process",
                "json_check_process", "check_analytic_results_process", "contact_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format (or no processes are specified)
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        info_msg  = "Using the old way to create the processes, this was removed!\n"
                        info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                        info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                        info_msg += "for a description of the new format"
                        raise Exception("StructuralMechanicsAnalysis: " + info_msg)

            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        raise Exception("Mixing of process initialization is not allowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                info_msg  = "Using the old way to create the gid-output, this was removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                raise Exception("StructuralMechanicsAnalysis: " + info_msg)
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _GetSimulationName(self):
        return "::[KSM Simulation]:: "

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
