from importlib import import_module
from typing import Any

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(models: 'dict[Kratos.Model]', parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"KratosCoSimAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"KratosCoSimAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return KratosCoSimAnalysisExecutionPolicy(parameters["name"].GetString(), models, parameters["settings"], optimization_problem)

class KratosCoSimAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, models: 'dict[Kratos.Model]', parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "model_names" : [],
            "analysis_module"  : "KratosMultiphysics",
            "analysis_type"    : "",
            "analysis_settings": {},
            "surrogate_io_vars": {}
        }""")
        self.models = models
        self.parameters = parameters
        self.optimization_problem = optimization_problem
        self.parameters.ValidateAndAssignDefaults(default_settings)

        analysis_module = parameters["analysis_module"].GetString()
        self.analysis_type = parameters["analysis_type"].GetString()
        self.analysis_settings = parameters["analysis_settings"]

        if analysis_module == "KratosMultiphysics":
            self.analysis_full_module, self.analysis_type = GetClassModuleFromKratos(self.analysis_type)
        else:
            self.analysis_full_module = f"{analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(self.analysis_type)}"

        self.analysis: AnalysisStage = getattr(import_module(self.analysis_full_module), self.analysis_type)(self.analysis_settings.Clone(), self.models)

        self.surrogate_io_settings = self.parameters["surrogate_io_vars"]
        self.variable_utils = Kratos.VariableUtils()

    def GetAnalysisModelPart(self):
        return self.analysis._GetSolver().GetComputingModelPart()

    def Initialize(self):
        self.analysis.Initialize()
        self.surrogate_io_initialize()

        # initialize model parts
        #self.model_parts = [self.model[model_part_name] for model_part_name in self.parameters["model_part_names"].GetStringArray()]

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        self.analysis.Finalize()

    def Execute(self):
        for k,v in self.analysis.models.items():
            model_part_names = v.GetModelPartNames()
            for model_part in model_part_names:
                mdp = v.GetModelPart(model_part)
                mdp.ProcessInfo.SetValue(Kratos.IS_RESTARTED, False)
                mdp.ProcessInfo.SetValue(Kratos.STEP, 0)
                mdp.ProcessInfo.SetValue(Kratos.TIME, 0.0)
                mdp.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0.0)
                if k == 'fluid':
                    mdp.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0.1)
        self.analysis.time = 0.0
        self.analysis.step = 0
        self.analysis.RunSolutionLoop()

        #self.Output()

        #self._OutputAnalysisData()

    def surrogate_io_train_append(self,control_var_val, step):
        #control = self.surrogate_io_settings["controlVars"].GetStringArray()
        sol_vars_info = self.surrogate_io_settings["solutionVars"].GetStringArray()
        nextSample = [step]
        #for element in control_vars_info:
        #val = eval("self.analysis_settings" + control_vars_info[1] + ".GetDouble()")
        nextSample.append(control_var_val)
        for element in sol_vars_info:
            element = eval(element)
            modelpart = self.models[element[2]][element[1]]
            if type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.DoubleVariable:
                var = Kratos.KratosGlobals.GetVariable(element[0])
                for node in modelpart.Nodes:
                    nextSample.append(node.GetSolutionStepValue(var))
            elif type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.Array1DVariable3:
                varX = Kratos.KratosGlobals.GetVariable(element[0] + '_X')
                varY = Kratos.KratosGlobals.GetVariable(element[0] + '_Y')
                for node in modelpart.Nodes:
                    valX = node.GetSolutionStepValue(varX)
                    valY = node.GetSolutionStepValue(varY)
                    nextSample.append(valX)
                    nextSample.append(valY)
        self.surrogate_io_train.append(nextSample)

    def surrogate_io_valid_append(self,control_var_val, step):
        #control_vars_info = self.surrogate_io_settings["controlVars"].GetStringArray()
        sol_vars_info = self.surrogate_io_settings["solutionVars"].GetStringArray()
        nextSample = [step]
        #for element in control_vars_info:
        #val = eval("self.analysis_settings" + control_vars_info[1] + ".GetDouble()")
        nextSample.append(control_var_val)
        for element in sol_vars_info:
            element = eval(element)
            modelpart = self.models[element[2]][element[1]]
            if type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.DoubleVariable:
                var = Kratos.KratosGlobals.GetVariable(element[0])
                for node in modelpart.Nodes:
                    nextSample.append(node.GetSolutionStepValue(var))
            elif type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.Array1DVariable3:
                varX = Kratos.KratosGlobals.GetVariable(element[0] + '_X')
                varY = Kratos.KratosGlobals.GetVariable(element[0] + '_Y')
                for node in modelpart.Nodes:
                    valX = node.GetSolutionStepValue(varX)
                    valY = node.GetSolutionStepValue(varY)
                    nextSample.append(valX)
                    nextSample.append(valY)
        self.surrogate_io_valid.append(nextSample)
        #self.surrogate_io_valid.append(header)
#
        #with open('readme.txt', 'w') as f:
#
        #    for val in values:
#
        #        if type(variable) == Kratos.Array1DVariable3:
        #            out = " " + " ".join( "{:e}".format(v) for v in val ) + "\n"
#
        #        f.write(out)
        #f.close()
        #a = 1

    def surrogate_io_initialize(self):
        self.surrogate_io_train = []
        self.surrogate_io_valid = []
        control_vars_info = self.surrogate_io_settings["controlVars"].GetStringArray()
        sol_vars_info = self.surrogate_io_settings["solutionVars"].GetStringArray()
        header = ["DoeNo"]
        #for element in control_vars_info:
        element = eval(control_vars_info[0])
        header.append(element)
        for element in sol_vars_info:
            element = eval(element)
            modelpart = self.models[element[2]][element[1]]
            if type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.DoubleVariable:
                for node in modelpart.Nodes:
                    header.append(element[0] + '_' + str(node.Id))
            elif type(Kratos.KratosGlobals.GetVariable(element[0])) == Kratos.Array1DVariable3:
                for node in modelpart.Nodes:
                    header.append(element[0] + '_X_' + str(node.Id))
                    header.append(element[0] + '_Y_' + str(node.Id))
        self.surrogate_io_train.append(header)
        self.surrogate_io_valid.append(header)

    def _OutputAnalysisData(self):
        unbuffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        for model_part in self.model_parts:
            for variable in self.nodal_solution_step_data_variables:
                nodal_field = Kratos.Expression.NodalExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(nodal_field, variable, True)
                unbuffered_data.SetValue(variable.Name(), nodal_field.Clone(), overwrite=True)
            for variable in self.nodal_data_value_variables:
                nodal_field = Kratos.Expression.NodalExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(nodal_field, variable, False)
                unbuffered_data.SetValue(variable.Name(), nodal_field.Clone(), overwrite=True)
            for variable in self.element_data_value_variables:
                elem_field = Kratos.Expression.ElementExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(elem_field, variable)
                unbuffered_data.SetValue(variable.Name(), elem_field.Clone(), overwrite=True)
            for variable in self.condition_data_value_variables:
                cond_field = Kratos.Expression.ConditionExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(cond_field, variable)
                unbuffered_data.SetValue(variable.Name(), cond_field.Clone(), overwrite=True)

    @staticmethod
    def __GetVariablesList(variable_names_list: 'list[str]') -> 'list[Any]':
        return [Kratos.KratosGlobals.GetVariable(variable_name) for variable_name in variable_names_list]



