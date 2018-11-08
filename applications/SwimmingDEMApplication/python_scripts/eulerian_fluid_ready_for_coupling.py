from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7


import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication

from fluid_dynamics_analysis import FluidDynamicsAnalysis

class python_parameters:
    def __init__(self):
        pass

class DEMCoupledFluidDynamicsAnalysis(FluidDynamicsAnalysis):

    def __init__(self, model, parameters=None):

        self.model = model
        self.pp = python_parameters()
        with open("ProjectParameters.json",'r') as parameter_file:
            self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
            gid_output_options = self.parameters["output_processes"]["gid_output"][0]["Parameters"]
            result_file_configuration = gid_output_options["postprocess_parameters"]["result_file_configuration"]
            nodal_results = result_file_configuration["nodal_results"]
            gauss_point_results = result_file_configuration["gauss_point_results"]
            self.pp.nodal_results = [nodal_results[i].GetString() for i in range(nodal_results.size())]
            self.pp.gauss_points_results = [gauss_point_results[i].GetString() for i in range(gauss_point_results.size())]
        self.pp.fluid_parameters = self.parameters

        super(DEMCoupledFluidDynamicsAnalysis, self).__init__(model, self.parameters)
        self.fluid_model_part = self._GetSolver().main_model_part

    def Initialize(self):
        self.AddFluidVariablesBySwimmingDEMAlgorithm()
        super(DEMCoupledFluidDynamicsAnalysis, self).Initialize()

    def AddFluidVariablesBySwimmingDEMAlgorithm(self):
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.pp.fluid_vars)

    def RunSingleTimeStep(self):
        self.InitializeSolutionStep()
        self._GetSolver().Predict()
        self._GetSolver().SolveSolutionStep()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

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
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = DEMCoupledFluidDynamicsAnalysis(model,parameters)
    simulation.Run()
