from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SM
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSM

# Other imports
import sys

# Import the base structural analysis
from structural_mechanics_analysis import StructuralMechanicsAnalysis as BaseClass

class ContactStructuralMechanicsAnalysis(BaseClass):
    """
    This class is the main-script of the ContactStructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Construct the base analysis.
        super(ContactStructuralMechanicsAnalysis, self).__init__(model, project_parameters)

    def Initialize(self):
        """ Initializing the Analysis """
        super(ContactStructuralMechanicsAnalysis, self).Initialize()
        self._GetSolver().SetEchoLevel(self.echo_level)
        # To avoid many prints
        if self.echo_level == 0:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """

        # First we check if one of the output processes will print output in this step this is done to save computation in case none of them will print
        is_output_step = False
        for output_process in self._GetListOfOutputProcesses():
            if output_process.IsOutputStep():
                is_output_step = True
                break

        if is_output_step:
            # Informing the output will be created
            KM.Logger.PrintWarning(self.__get_simulation_name_output(), "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KM.STEP])
            KM.Logger.PrintWarning(self.__get_simulation_name_output(), "TIME: ", self.time)

        # Creating output
        super(ContactStructuralMechanicsAnalysis, self).OutputSolutionStep()

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        import python_solvers_wrapper_contact_structural
        return python_solvers_wrapper_contact_structural.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        return "::[KCSM Simulation]:: "

    def __get_simulation_name_output(self):
        return "::[KCSM Simulation: Output]:: "

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 contact_structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ContactStructuralMechanicsAnalysis(model, parameters)
    simulation.Run()
