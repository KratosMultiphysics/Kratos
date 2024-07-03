from KratosMultiphysics import Model, Parameters
import KratosMultiphysics.FluidDynamicsApplication

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis_rve import FluidDynamicsAnalysisRve

class DEMCoupledFluidDynamicsAnalysisPeriodicity(FluidDynamicsAnalysisRve):

    def __init__(self, model, parameters=None, variables_management=None):
        self.model = model
        self.sdem_project_parameters = parameters
        self.project_parameters = self.sdem_project_parameters['fluid_parameters']
        self.vars_man = variables_management

        super().__init__(model, self.project_parameters)
        self.fluid_model_part = self._GetSolver().main_model_part

    def Initialize(self):
        self.AddFluidVariablesForSwimmingDEM()
        super().Initialize()

    def AddFluidVariablesForSwimmingDEM(self):
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.vars_man.fluid_vars)

    # Ensures that the avg velocity in the RVE domain remains 0: in this case, we do not want to, if not, call the base method
    def _AvgVelocityCorrection(self):
        pass

    def _CreateSolver(self):
        if (self.project_parameters["solver_settings"]["solver_type"].GetString() == "MonolithicDEM"):
            from KratosMultiphysics.SwimmingDEMApplication.python_solvers_wrapper_fluidDEM import CreateSolver
            return CreateSolver(self.model, self.project_parameters)
        else:
            from KratosMultiphysics.FluidDynamicsApplication.python_solvers_wrapper_fluid import CreateSolver
            return CreateSolver(self.model, self.project_parameters)

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
