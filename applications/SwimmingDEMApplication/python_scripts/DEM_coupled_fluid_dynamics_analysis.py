from KratosMultiphysics import Model, Parameters
import KratosMultiphysics.FluidDynamicsApplication

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class DEMCoupledFluidDynamicsAnalysis(FluidDynamicsAnalysis):

    def __init__(self, model, parameters=None, variables_management=None):
        self.model = model
        self.sdem_project_parameters = parameters
        self.project_parameters = self.sdem_project_parameters['fluid_parameters']
        self.vars_man = variables_management

        super(DEMCoupledFluidDynamicsAnalysis, self).__init__(model, self.project_parameters)
        self.fluid_model_part = self._GetSolver().main_model_part

    def Initialize(self):
        self.AddFluidVariablesForSwimmingDEM()
        super(DEMCoupledFluidDynamicsAnalysis, self).Initialize()

    def AddFluidVariablesForSwimmingDEM(self):
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.vars_man.fluid_vars)

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Parameters(parameter_file.read())

    model = Model()
    simulation = DEMCoupledFluidDynamicsAnalysis(model, parameters)
    simulation.Run()
