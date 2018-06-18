from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.PoromechanicsApplication

from poromechanics_analysis import PoromechanicsAnalysis

#TODO

class PoromechanicsFractureAnalysis(PoromechanicsAnalysis):
    '''Main script for poromechanics simulations.'''

    def __init__(self,model,parameters):
        super(PoromechanicsFractureAnalysis,self).__init__(model,parameters)

    def Initialize(self):
        '''
        Construct and initialize all classes and tools used in the simulation loop.
        '''

        super(PoromechanicsFractureAnalysis,self).Initialize()

        # Initialize Fracture Propagation Utility
        import poromechanics_fracture_propagation_utility
        self.fracture_utility = poromechanics_fracture_propagation_utility.FracturePropagationUtility(self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE],
                                                                                                    self.project_parameters["problem_data"]["problem_name"].GetString(),
                                                                                                    self.project_parameters["solver_settings"]["move_mesh_flag"].GetBool())

    def OutputSolutionStep(self):

        super(PoromechanicsFractureAnalysis,self).OutputSolutionStep()

        # Check Fracture Propagation Utility
        if self.fracture_utility.IsPropagationStep():
            self.model,self.main_model_part,self.solver,self.list_of_processes,self.output = self.fracture_utility.CheckPropagation(self.model,
                                                                                                                        self.main_model_part,
                                                                                                                        self.solver,
                                                                                                                        self.list_of_processes,
                                                                                                                        self.output)

    def Finalize(self):

        # Finalize Fracture Propagation Utility
        self.fracture_utility.Finalize()

        super(PoromechanicsFractureAnalysis,self).Finalize()

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python poromechanics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python poromechanics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = PoromechanicsFractureAnalysis(model,parameters)
    simulation.Run()
