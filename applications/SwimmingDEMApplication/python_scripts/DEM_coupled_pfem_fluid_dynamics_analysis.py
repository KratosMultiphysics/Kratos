from KratosMultiphysics import Model, Parameters
import KratosMultiphysics.PfemFluidDynamicsApplication
from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis import PfemFluidDynamicsAnalysis

from importlib import import_module

class DEMCoupledPFEMFluidDynamicsAnalysis(PfemFluidDynamicsAnalysis):

    def __init__(self, model, parameters=None, variables_management=None):
        self.model = model
        self.sdem_project_parameters = parameters
        self.project_parameters = self.sdem_project_parameters['fluid_parameters']
        self.dimension = self.project_parameters["solver_settings"]["domain_size"].GetInt()
        self.vars_man = variables_management
        variables_management.nodal_results, variables_management.gauss_points_results = [], []

        if self.project_parameters.Has('output_configuration'):
            gid_output_options = self.project_parameters["output_configuration"]["result_file_configuration"]
            gauss_point_results = gid_output_options["gauss_point_results"]
            nodal_variables = self.project_parameters["output_configuration"]["result_file_configuration"]["nodal_results"]
            variables_management.nodal_results = [nodal_variables[i].GetString() for i in range(nodal_variables.size())]
            variables_management.gauss_points_results = [gauss_point_results[i].GetString() for i in range(gauss_point_results.size())]

        super().__init__(model, self.project_parameters)
        self.fluid_model_part = self._GetSolver().main_model_part

    def Initialize(self):
        self.AddFluidVariablesForSwimmingDEM()
        super().Initialize()

    def AddFluidVariablesForSwimmingDEM(self):
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.vars_man.fluid_vars)

    def RunSingleTimeStep(self):
        self.InitializeSolutionStep()
        self._GetSolver().Predict()
        self._GetSolver().SolveSolutionStep()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

    def _CreateSolver(self):
        python_module_name = "KratosMultiphysics.SwimmingDEMApplication"
        full_module_name = python_module_name + "." + self.project_parameters["solver_settings"]["solver_type"].GetString()
        solver_module = import_module(full_module_name)
        solver = solver_module.CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver

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
    simulation = DEMCoupledFluidDynamicsAnalysis(model,parameters)
    simulation.Run()