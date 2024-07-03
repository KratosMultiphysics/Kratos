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
        # self._CreateModelers()
        # self._ModelersSetupGeometryModel()
        # self._ModelersSetupGeometryModel()
        # self._ModelersSetupModelPart()

        # self._GetSolver().ImportModelPart()
        # self._GetSolver().PrepareModelPart()
        # self._GetSolver().AddDofs()

        # self.ModifyInitialProperties()
        # self.ModifyInitialGeometry()

        # ##here we initialize user-provided processes
        # self.__CreateListOfProcesses() # has to be done after importing and preparing the ModelPart
        # for process in self._GetListOfProcesses():
        #     process.ExecuteInitialize()

        # self._GetSolver().Initialize()

        # self.Check()

        # self.ModifyAfterSolverInitialize()

        # for process in self._GetListOfProcesses():
        #     process.ExecuteBeforeSolutionLoop()

        # ## Stepping and time settings
        # self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        # if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
        #     self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        # else:
        #     self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
        #     self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time

        # ## If the echo level is high enough, print the complete list of settings used to run the simulation
        # if self.echo_level > 1:
        #     with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
        #         parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        # KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def AddFluidVariablesForSwimmingDEM(self):
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.vars_man.fluid_vars)

    def _CreateSolver(self):
        if (self.project_parameters["solver_settings"]["solver_type"].GetString() == "MonolithicDEM"):
            from KratosMultiphysics.SwimmingDEMApplication.python_solvers_wrapper_fluidDEM import CreateSolver
            return CreateSolver(self.model, self.project_parameters)
        else:
            from KratosMultiphysics.FluidDynamicsApplication.python_solvers_wrapper_fluid import CreateSolver
            return CreateSolver(self.model, self.project_parameters)

    def __CreateListOfProcesses(self):
        """This function creates the processes and the output-processes
        """
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)
        deprecated_output_processes    = self._CheckDeprecatedOutputProcesses(self._list_of_processes)
        order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
        self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
        self._list_of_processes.extend(self._list_of_output_processes) # Adding the output processes to the regular processes
        self._list_of_output_processes.extend(deprecated_output_processes)

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
