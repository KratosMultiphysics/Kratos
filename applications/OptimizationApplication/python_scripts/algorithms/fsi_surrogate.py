import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.opt_convergence import CreateConvergenceCriteria
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAlgorithmTimeLogger
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
import pandas as pd
import numpy as np

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return FsiSurrogate(model, parameters, optimization_problem)

class FsiSurrogate(Algorithm):
    """
        A classical steepest descent algorithm to solve unconstrainted optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "controls"          : [],
            "echo_level"        : 0,
            "surrogate_settings": {}
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        
        self.master_control = MasterControl()
        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)

        #settings = parameters["settings"]
        #settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.echo_level = 0

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        pass

    @time_decorator()
    def Initialize(self):
        self.converged = False
        self.step = 0
        self._import_DOE()
        self._assign_doe_to_controls()

    def Finalize(self):
        self.__objective.Finalize()
        self.master_control.Finalize()

    @time_decorator()
    def UpdateControlTrain(self) -> None:
        step = self.step
        for control in self.master_control.GetListOfControls():
            control.NextTrainingSet(step)

    def UpdateControlValid(self) -> None:
        step = self.step
        for control in self.master_control.GetListOfControls():
            control.NextValidationSet(step)

    @time_decorator()
    def OutputTrainingData(self):
        for output_process in self._optimization_problem.GetListOfProcesses("output_processes"):
            output_process.PrintOutputTrainingSamples(self.step)
    
    def OutputValidationData(self):
        for output_process in self._optimization_problem.GetListOfProcesses("output_processes"):
            output_process.PrintOutputValidationSamples(self.step)

    @time_decorator()
    def Solve(self):
        nTrain = 0
        nValid = 0
        if self.trainDOE:
            nTrain = self.trainDOE["nTrain"]
        if self.validDOE:
            nValid = self.validDOE["nValid"]

        while self.step <= nTrain:
                self.UpdateControlTrain()
                self.fsisurrogate = self._optimization_problem.GetExecutionPolicy("fsi")
                self.fsisurrogate.Execute()
                self.fsisurrogate.surrogate_io_train_append(self.trainDOE['YOUNG_MODULUS'][self.step],self.step)

                self.OutputTrainingData()
#
                #self.converged = self.__convergence_criteria.IsConverged()
#
                self.AdvanceStep()

        self.ResetStep()

        while self.step <= nValid:
            self.UpdateControlValid()
            self.fsisurrogate = self._optimization_problem.GetExecutionPolicy("fsi")
            self.fsisurrogate.Execute()
            self.fsisurrogate.surrogate_io_valid_append(self.validDOE['YOUNG_MODULUS'][self.step],self.step)
            self.OutputValidationData()
            self.AdvanceStep()

        self.print_excel()
        
        return self.converged
    
    def _import_DOE(self):
        if self.parameters["surrogate_settings"]["DOE_csv_files"].Has("train"):
            self.trainDOE = pd.read_csv(self.parameters["surrogate_settings"]["DOE_csv_files"]["train"].GetString()).to_dict('list')
            listVals = list(self.trainDOE.values())
            self.trainDOE["nTrain"] = len(listVals[0]) - 1
            print("Training DOE file read complete",'\n')
        if self.parameters["surrogate_settings"]["DOE_csv_files"].Has("valid"):
            self.validDOE = pd.read_csv(self.parameters["surrogate_settings"]["DOE_csv_files"]["valid"].GetString()).to_dict('list')
            listVals = list(self.validDOE.values())
            self.validDOE["nValid"] = len(listVals[0]) - 1
            print("Validation DOE file read complete",'\n')

    def _assign_doe_to_controls(self):
        listOfControls = self.master_control.GetListOfControls()
        for control in listOfControls:
            variable = control.GetControlVarName()
            if self.trainDOE:
                values = self.trainDOE[variable]
                control.SetTrainDOEValues(values)
            if self.validDOE:
                values = self.validDOE[variable]
                control.SetValidDOEValues(values)
    
    def AdvanceStep(self):
        current_step = self.step
        self.step = current_step + 1

    def ResetStep(self):
        self.step = 0
    
    def print_excel(self):
        io_train_df = pd.DataFrame(self.fsisurrogate.surrogate_io_train)
        io_train_df.to_excel(self.parameters["surrogate_settings"]["surrogate_io_train_excel"].GetString(),header=False, index=False)

        io_valid_df = pd.DataFrame(self.fsisurrogate.surrogate_io_valid)
        io_valid_df.to_excel(self.parameters["surrogate_settings"]["surrogate_io_valid_excel"].GetString(),header=False, index=False)
