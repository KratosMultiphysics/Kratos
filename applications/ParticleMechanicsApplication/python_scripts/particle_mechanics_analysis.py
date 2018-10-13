from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as ParticleMechanicsApplication

# Importing the solvers (if available)
try:
    import KratosMultiphysics.ExternalSolversApplication
    KratosMultiphysics.Logger.PrintInfo("ExternalSolversApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("ExternalSolversApplication", "not imported")

# Importing the base class
from analysis_stage import AnalysisStage

class ParticleMechanicsAnalysis(AnalysisStage):
    """
    This class is the main-script of the ParticleMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initializing the parameters
        solver_settings = project_parameters["solver_settings"]

        if solver_settings.Has("domain_size") and project_parameters["problem_data"].Has("domain_size"):
            warn_msg  = '"domain_size" defined both in "problem_data" and "solver_settings"!'
            warn_msg += 'the definition in the "solver_settings" will be employed.'
            KratosMultiphysics.Logger.PrintWarning("ParticleMechanicsAnalysis", warn_msg)

        if solver_settings.Has("model_part_name") and project_parameters["problem_data"].Has("model_part_name"):
            warn_msg  = '"model_part_name" defined both in problem_data" and "solver_settings"!'
            warn_msg += 'the definition in the "solver_sett"ings" will be employed.'
            KratosMultiphysics.Logger.PrintWarning("ParticleMechanicsAnalysis", warn_msg)

        if solver_settings.Has("time_stepping") and project_parameters["problem_data"].Has("time_Step"):
            warn_msg  = '"time_stepping" defined both in "problem_data" and "solver_settings"!'
            warn_msg += 'the definition in the "solver_settings" will be employed.'
            KratosMultiphysics.Logger.PrintWarning("ParticleMechanicsAnalysis", warn_msg)

        if not solver_settings.Has("time_stepping"):
            KratosMultiphysics.Logger.PrintWarning("ParticleMechanicsAnalysis", "Using the old way to pass the time_step, this will be removed!")
            time_stepping_params = KratosMultiphysics.Parameters("{}")
            time_stepping_params.AddValue("time_step", project_parameters["problem_data"]["time_step"])
            solver_settings.AddValue("time_stepping", time_stepping_params)

        if not solver_settings.Has("domain_size"):
            KratosMultiphysics.Logger.PrintWarning("ParticleMechanicsAnalysis", "Using the old way to pass the domain_size, this will be removed!")
            solver_settings.AddEmptyValue("domain_size")
            solver_settings["domain_size"].SetInt(project_parameters["problem_data"]["domain_size"].GetInt())

        if not solver_settings.Has("model_part_name"):
            KratosMultiphysics.Logger.PrintWarning("ParticleMechanicsAnalysis", "Using the old way to pass the model_part_name, this will be removed!")
            solver_settings.AddEmptyValue("model_part_name")
            solver_settings["model_part_name"].SetString(project_parameters["problem_data"]["model_part_name"].GetString())

        # Import parallel modules if needed
        # has to be done before the base-class constuctor is called (in which the solver is constructed)
        # TODO: currently MPI parallelization is not present in ParticleMechanicsApplication
        # TODO: the following import lines will be kept here for future reference
        if (project_parameters["problem_data"]["parallel_type"].GetString() == "MPI"):
            warn_msg  = 'Currently MPI parallelization is not present in ParticleMechanicsApplication!'
            KratosMultiphysics.Logger.PrintWarning("ParticleMechanicsAnalysis", warn_msg)
            # import KratosMultiphysics.MetisApplication as MetisApplication
            # import KratosMultiphysics.TrilinosApplication as TrilinosApplication

        super(ParticleMechanicsAnalysis, self).__init__(model, project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        import python_solvers_wrapper_particle
        return python_solvers_wrapper_particle.CreateSolver(self.model, self.project_parameters)
    
    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes"""
        list_of_processes = super(ParticleMechanicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["constraints_process_list", "loads_process_list", "list_other_processes", "gravity"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                info_msg  = "Using the old way to create the processes, this will be removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                KratosMultiphysics.Logger.PrintWarning("ParticleMechanicsAnalysis", info_msg)
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
                KratosMultiphysics.Logger.PrintInfo("ParticleMechanicsAnalysis", info_msg)
                gid_output= self._SetUpGiDOutput()
                list_of_processes += [gid_output,]
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _SetUpGiDOutput(self):
        '''Initialize a GiD output instance'''
        


if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 particle_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 particle_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ParticleMechanicsAnalysis(model, parameters)
    simulation.Run()
