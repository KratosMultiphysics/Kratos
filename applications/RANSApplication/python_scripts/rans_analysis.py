from sys import argv

import KratosMultiphysics as Kratos
from KratosMultiphysics import IsDistributedRun

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class RANSAnalysis(FluidDynamicsAnalysis):
    '''Main script for fluid dynamics simulations using the navier_stokes family of python solvers.'''

    def __init__(self,model,parameters):
        super().__init__(model, parameters)
        self.delta_time = 0.0

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """

        while self.KeepAdvancingSolutionLoop():
            new_time = self._GetSolver().AdvanceInTime(self.time)
            self.delta_time = new_time - self.time
            self.time = new_time
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def _ReInitializeSolver(self):
        """ This reinitializes after remesh """
        self._GetSolver().Clear()

        self._GetSolver().PrepareModelPart()
        self._GetSolver().AddDofs()

        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        ##here we initialize user-provided processes
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()
        self.Check()

        self.ModifyAfterSolverInitialize()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Reinitialized solver after re-mesh")

    def _CreateSolver(self):
        solver_type = self.project_parameters["solver_settings"]["solver_type"].GetString()
        if IsDistributedRun():
            parallelism = "MPI"
        else:
            parallelism = "OpenMP"

        if (solver_type == "CoupledRANS"):
            from KratosMultiphysics.RANSApplication.coupled_rans_solver import CoupledRANSSolver
            return CoupledRANSSolver(self.model, self.project_parameters["solver_settings"])
        elif (solver_type == "CoupledRANSALE"):
            from KratosMultiphysics.RANSApplication.coupled_rans_ale_solver import CoupledRANSALESolver
            return CoupledRANSALESolver(self.model, self.project_parameters["solver_settings"], parallelism)
        else:
            raise Exception("Unknown solver type requested. [ solver_type = " + solver_type + " ].")

    def _GetSimulationName(self):
        return "RANS Analysis"

    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        return self.end_time - self.time > 0.5 * self.delta_time

if __name__ == '__main__':
    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python rans_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python rans_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = RANSAnalysis(model,parameters)
    simulation.Run()
