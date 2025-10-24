import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM
import csv
import math
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.SystemIdentificationApplication as KratosDT
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.response_sensitivity_analysis import ResponseSensitivityAnalysis
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.sensor_sensitivity_adjoint_dynamic_solver import SensorSensitivityAdjointDynamicSolver
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _):
    return SystemIdentificationDynamicAnalysis(model, parameters["settings"])

class SystemIdentificationDynamicAnalysis(ResponseSensitivityAnalysis):
    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

        # compute delta time
        delta_time = self._GetSolver()._ComputeDeltaTime()
        if delta_time >= 0.0:
            raise Exception("Adjoints are solved in reverse in time. Hence, delta_time should be negative. [ delta_time = " + str(delta_time) + "].")

        # update start time and end time. Adjoints are solved reverse in time,
        # hence we include an additional time step at the end
        # to properly initialize initial conditions for adjoint problem
        # and start time is also shifted by one to keep the same number
        # of steps
        delta_time *= -1.0
        problem_data = project_parameters["problem_data"]
        self.start_time = problem_data["start_time"].GetDouble() + delta_time
        self.end_time = problem_data["end_time"].GetDouble() + delta_time

        # now set start_time, end_time of the problem_data bt flipping them
        self.project_parameters["problem_data"]["start_time"].SetDouble(self.end_time)
        self.project_parameters["problem_data"]["end_time"].SetDouble(self.start_time)  

        self.simulation_time = self.end_time
        print("simul time ", self.simulation_time)    
        self.num_of_timesteps = math.ceil((self.end_time - self.start_time) / delta_time)
        print("num of measurements are ", self.num_of_timesteps)
        if self.num_of_timesteps <= 0:
            raise RuntimeError(f"Number of measurements are {self.num_of_timesteps}. Make sure they are > 0.")

    def NumberOfTimeSteps(self) -> int:
        return self.num_of_timesteps

    def _CreateSolver(self) -> SensorSensitivityAdjointDynamicSolver:
        return SensorSensitivityAdjointDynamicSolver(self.model, self.project_parameters["solver_settings"])

    def _GetSimulationName(self) -> str:
        return "::[SystemIdentificationAnalysis]:: "

    def CalculateGradient(self, response_function: Kratos.AdjointResponseFunction) -> 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]':
        optimization_step = self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP]

        Kratos.VariableUtils().ClearHistoricalData(self._GetSolver().GetComputingModelPart())
        
        print("In CalcGrad, Step (Opt step is) ", optimization_step)
        self._GetSolver().SetResponseFunction(response_function)
 
        self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP] = self.num_of_timesteps + 1
        self.time =  self.simulation_time
        self._GetSolver().main_model_part.CloneTimeStep(self.project_parameters["problem_data"]["start_time"].GetDouble())
 
        while self.KeepAdvancingSolutionLoop():

            self.time = self._AdvanceTime()
            print("loop time is ", self.time)
            self.InitializeSolutionStep()
            response_function.InitializeSolutionStep()
            self._GetSolver().Predict()

            node1 : Kratos.Node = self._GetSolver().GetComputingModelPart().GetNode(1)
            
            if not node1.GetDof(KratosSM.ADJOINT_DISPLACEMENT_X).IsFixed():
                raise RuntimeError("Not fixed node 1")

            is_converged = self._GetSolver().SolveSolutionStep()
            self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            response_function.FinalizeSolutionStep()
            self.FinalizeSolutionStep()
            
            self.OutputSolutionStep()

            
        self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP] = optimization_step
        gradients = self._GetSolver().GetSensitivities()
        return gradients
    
    def InitializeSolutionStep(self):

        exec_policy = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosDT.TEST_ANALYSIS_NAME]
        step = self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP]

        super().InitializeSolutionStep()
        
     
    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        Note that as the adjoint problem is solved backward in time, the stopping criteria is current time being larger than the start one
        """
        print("Inside Advancining loop check step ", self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP])
        return self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP] > 1

    def PrintAnalysisStageProgressInformation(self):
        process_info = self._GetSolver().GetComputingModelPart().ProcessInfo
        Kratos.Logger.PrintInfo(self._GetSimulationName(), f"Step {process_info[Kratos.STEP]}: Computed sensitivities for response \"{process_info[KratosSI.SENSOR_NAME]}\" using \"{process_info[KratosSI.TEST_ANALYSIS_NAME]}\" analysis.")

    def __CheckIfSolveSolutionStepReturnsAValue(self, is_converged):
        """In case the solver does not return the state of convergence
        (same as the SolvingStrategy does) then issue ONCE a deprecation-warning

        """
        if is_converged is None:
            if not hasattr(self, '_map_ret_val_depr_warnings'):
                self._map_ret_val_depr_warnings = []
            solver_class_name = self._GetSolver().__class__.__name__
            # used to only print the deprecation-warning once
            if not solver_class_name in self._map_ret_val_depr_warnings:
                self._map_ret_val_depr_warnings.append(solver_class_name)
                warn_msg  = 'Solver "{}" does not return '.format(solver_class_name)
                warn_msg += 'the state of convergence from "SolveSolutionStep"'
                IssueDeprecationWarning("AnalysisStage", warn_msg)

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = SystemIdentificationDynamicAnalysis(model, parameters)
    simulation.Run()
