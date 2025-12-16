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

def Factory(models: 'dict[Kratos.Model]', parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return FsiSurrogate(models, parameters, optimization_problem)

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

    def __init__(self, models: 'dict[Kratos.Model]', parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.models = models
        self.parameters = parameters
        self.surrogate_settings = parameters["surrogate_settings"]
        self._optimization_problem = optimization_problem

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        
        self.master_control = MasterControl()
        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)
    
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
        self._surrogate_io_initialize()

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
    def OutputTrainingData(self,f):
        nextSample = self.surrogate_io_train_append()
        dfNextSample = pd.DataFrame([nextSample])
        dfNextSample.to_csv(f,index=False,header=False,mode='a')
        f.flush()
        for output_process in self._optimization_problem.GetListOfProcesses("output_processes"):
            output_process.PrintOutputTrainingSamples(self.step)
    
    def OutputValidationData(self,f):
        nextSample = self.surrogate_io_valid_append()
        dfNextSample = pd.DataFrame([nextSample])
        dfNextSample.to_csv(f,index=False,header=False,mode='a')
        f.flush()
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

        with open(self.surrogate_settings["samples_csv"]["train"].GetString(), 'a') as trainSamplesCSV:
            while self.step <= nTrain:
                    self.UpdateControlTrain()
                    self.fsisurrogate = self._optimization_problem.GetExecutionPolicy("fsi")
                    self.fsisurrogate.Execute()
                    self.OutputTrainingData(trainSamplesCSV)
                    self.AdvanceStep()
            trainSamplesCSV.close()

        self.ResetStep()

        with open(self.surrogate_settings["samples_csv"]["valid"].GetString(), 'a') as validSamplesCSV:
            while self.step <= nValid:
                self.UpdateControlValid()
                self.fsisurrogate = self._optimization_problem.GetExecutionPolicy("fsi")
                self.fsisurrogate.Execute()
                self.OutputValidationData(validSamplesCSV)
                self.AdvanceStep()
            validSamplesCSV.close()
        
        return self.converged
    
    def _import_DOE(self):
        if self.parameters["surrogate_settings"]["DOE_csv"].Has("train"):
            self.trainDOE = pd.read_csv(self.parameters["surrogate_settings"]["DOE_csv"]["train"].GetString()).to_dict('list')
            listVals = list(self.trainDOE.values())
            self.trainDOE["nTrain"] = len(listVals[0]) - 1
            print("Training DOE file read complete",'\n')
        if self.parameters["surrogate_settings"]["DOE_csv"].Has("valid"):
            self.validDOE = pd.read_csv(self.parameters["surrogate_settings"]["DOE_csv"]["valid"].GetString()).to_dict('list')
            listVals = list(self.validDOE.values())
            self.validDOE["nValid"] = len(listVals[0]) - 1
            print("Validation DOE file read complete",'\n')

    def _assign_doe_to_controls(self):
        listOfControls = self.master_control.GetListOfControls()
        for control in listOfControls:
            variables = control.GetControlVarNames()
            if self.trainDOE:
                doe = { k: v for (k,v) in self.trainDOE.items() if k in variables}
                control.SetTrainDOEValues(doe)
            if self.validDOE:
                doe = { k: v for (k,v) in self.validDOE.items() if k in variables}
                control.SetValidDOEValues(doe)
    
    def _surrogate_io_initialize(self):
        self.surrogate_io_train = []
        self.surrogate_io_valid = []
        control_vars_info = self.surrogate_settings["variables"]["controlVars"].GetStringArray()
        sol_vars_info = self.surrogate_settings["variables"]["solutionVars"].GetStringArray()
        header = ["DoeNo"]
        for element in control_vars_info:
            header.append(element)
        for element in sol_vars_info:
            element = eval(element)
            modelpart = self.models[element[2]][element[1]]
            if type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.DoubleVariable:
                if element[0] == "PRESSURE_COEFFICIENT" :
                    for cond in modelpart.Conditions:
                        header.append(element[0] + '_' + str(cond.Id))
                elif (element[0] == "LIFT_COEFFICIENT") or (element[0] == "DRAG_COEFFICIENT"):
                    header.append(element[0])
                else :
                    for node in modelpart.Nodes:
                        header.append(element[0] + '_' + str(node.Id))
            elif type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.Array1DVariable3:
                for node in modelpart.Nodes:
                    header.append(element[0] + '_X_' + str(node.Id))
                    header.append(element[0] + '_Y_' + str(node.Id))
                    header.append(element[0] + '_Z_' + str(node.Id))
        dfHeader = pd.DataFrame([header])
        dfHeader.to_csv(self.surrogate_settings["samples_csv"]["train"].GetString(),index=False, header=False)
        dfHeader.to_csv(self.surrogate_settings["samples_csv"]["valid"].GetString(),index=False, header=False)

    def surrogate_io_train_append(self):
        control_vars_info = self.surrogate_settings["variables"]["controlVars"].GetStringArray()
        sol_vars_info = self.surrogate_settings["variables"]["solutionVars"].GetStringArray()
        nextSample = [self.step]
        for element in control_vars_info:
            val = self.trainDOE[element][self.step] 
            nextSample.append(val)
        for element in sol_vars_info:
            element = eval(element)
            modelpart = self.models[element[2]][element[1]]
            root_model_part = modelpart.GetRootModelPart()
            if type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.DoubleVariable:
                var = Kratos.KratosGlobals.GetVariable(element[0])
                if element[0] == "PRESSURE_COEFFICIENT":
                      for cond in modelpart.Conditions:
                           nextSample.append(cond.GetValue(var))
                elif element[0] == "LIFT_COEFFICIENT" or element[0] == "DRAG_COEFFICIENT":
                    nextSample.append(root_model_part.ProcessInfo.GetValue(var))
                else:
                	for node in modelpart.Nodes:
                            nextSample.append(node.GetSolutionStepValue(var))
            elif type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.Array1DVariable3:
                varX = Kratos.KratosGlobals.GetVariable(element[0] + '_X')
                varY = Kratos.KratosGlobals.GetVariable(element[0] + '_Y')
                varZ = Kratos.KratosGlobals.GetVariable(element[0] + '_Z')
                if element[0] == "VELOCITY":
                    for node in modelpart.Nodes:
                        valX = node.GetValue(varX)
                        valY = node.GetValue(varY)
                        valZ = node.GetValue(varZ)
                        nextSample.append(valX)
                        nextSample.append(valY)
                        nextSample.append(valZ)
                else:
                    for node in modelpart.Nodes:
                        valX = node.GetSolutionStepValue(varX)
                        valY = node.GetSolutionStepValue(varY)
                        valZ = node.GetSolutionStepValue(varZ)
                        nextSample.append(valX)
                        nextSample.append(valY)
                        nextSample.append(valZ)
        return nextSample
    
    def surrogate_io_valid_append(self):
        control_vars_info = self.surrogate_settings["variables"]["controlVars"].GetStringArray()
        sol_vars_info = self.surrogate_settings["variables"]["solutionVars"].GetStringArray()
        nextSample = [self.step]
        for element in control_vars_info:
            val = self.validDOE[element][self.step] 
            nextSample.append(val)
        for element in sol_vars_info:
            element = eval(element)
            modelpart = self.models[element[2]][element[1]]
            root_model_part = modelpart.GetRootModelPart()
            if type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.DoubleVariable:
                var = Kratos.KratosGlobals.GetVariable(element[0])
                if element[0] == "PRESSURE_COEFFICIENT":
                    for cond in modelpart.Conditions:
                        nextSample.append(cond.GetValue(var))
                elif element[0] == "LIFT_COEFFICIENT" or element[0] == "DRAG_COEFFICIENT":
                    nextSample.append(root_model_part.ProcessInfo.GetValue(var))
                else:
                    for node in modelpart.Nodes:
                        nextSample.append(node.GetSolutionStepValue(var))
            elif type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.Array1DVariable3:
                varX = Kratos.KratosGlobals.GetVariable(element[0] + '_X')
                varY = Kratos.KratosGlobals.GetVariable(element[0] + '_Y')
                varZ = Kratos.KratosGlobals.GetVariable(element[0] + '_Z')
                if element[0] == "VELOCITY":
                    for node in modelpart.Nodes:
                        valX = node.GetValue(varX)
                        valY = node.GetValue(varY)
                        valZ = node.GetValue(varZ)
                        nextSample.append(valX)
                        nextSample.append(valY)
                        nextSample.append(valZ)
                else:
                    for node in modelpart.Nodes:
                        valX = node.GetSolutionStepValue(varX)
                        valY = node.GetSolutionStepValue(varY)
                        valZ = node.GetSolutionStepValue(varZ)
                        nextSample.append(valX)
                        nextSample.append(valY)
                        nextSample.append(valZ)
        return nextSample
    
    def AdvanceStep(self):
        current_step = self.step
        self.step = current_step + 1

    def ResetStep(self):
        self.step = 0
