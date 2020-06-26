from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

## Importing modules -----------------------------------------------------------------------------------------

# Import kratos core and applications
import KratosMultiphysics as Kratos

from KratosMultiphysics.FluidTransportApplication.fluid_transport_analysis import FluidTransportAnalysis
from KratosMultiphysics.analysis_stage import AnalysisStage

class FluidTransportCouplingAnalysis(FluidTransportAnalysis):
    '''Main script for fluid + fluid transport simulations.'''

    def __init__(self,model,parameters):

        # Creating solver and model part and adding variables
        super(FluidTransportCouplingAnalysis,self).__init__(model,parameters)

#TODO if we update processes input in fluid transport interface, this file will not be necessary

    def _CreateProcesses(self, parameter_name, initialization_order):
        return AnalysisStage._CreateProcesses(self,parameter_name,initialization_order)

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_transport_coupling_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_transport_coupling_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = FluidTransportCouplingAnalysis(model,parameters)
    simulation.Run()
