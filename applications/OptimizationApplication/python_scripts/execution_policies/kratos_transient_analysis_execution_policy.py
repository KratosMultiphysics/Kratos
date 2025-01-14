from importlib import import_module

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.analysis_stage import AnalysisStage

from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"KratosTransientAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"KratosTransientAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return KratosTransientAnalysisExecutionPolicy(parameters["name"].GetString(), model, parameters["settings"])

class KratosTransientAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        self.model = model

        default_settings = Kratos.Parameters("""{
            "model_part_names"  : [],
            "analysis_module"  : "KratosMultiphysics",
            "analysis_type"    : "",
            "analysis_settings": {}
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        self.analysis_module = parameters["analysis_module"].GetString()
        self.analysis_type = parameters["analysis_type"].GetString()
        self.analysis_settings = parameters["analysis_settings"]

        self.analysis_start_time = self.analysis_settings["problem_data"]["start_time"].GetDouble()
        self.analysis_end_time = self.analysis_settings["problem_data"]["end_time"].GetDouble()
        self.analysis_time_step_settings = self.analysis_settings["solver_settings"]["time_stepping"]

        if self.analysis_module == "KratosMultiphysics":
            self.analysis_full_module, self.analysis_type = GetClassModuleFromKratos(self.analysis_type)
        else:
            self.analysis_full_module = f"{self.analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(self.analysis_type)}"

        self.model_part_names = parameters["model_part_names"].GetStringArray()
        self.model_parts: 'list[Kratos.ModelPart]' = []
        self.analysis: AnalysisStage = getattr(import_module(self.analysis_full_module), self.analysis_type)(self.model, self.analysis_settings.Clone())
    
    def Initialize(self) -> None:
        self.analysis.Initialize()
        self.model_parts = [self.model[model_part_name] for model_part_name in self.model_part_names]
    
    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        # self.analysis.Finalize()
        pass

    def Execute(self):
        # Reset model part and analysis
        for model_part in self.model_parts:
            KratosOA.OptimizationUtils.ResetModelPartNodalSolutionStepData(model_part)
            # model_part.ProcessInfo[Kratos.TIME] = 0.0
            # model_part.ProcessInfo[Kratos.DELTA_TIME] = 0.0
        # self.analysis: AnalysisStage = getattr(import_module(self.analysis_full_module), self.analysis_type)(self.model, self.analysis_settings.Clone())
        self.GetAnalysisModelPart().ProcessInfo[Kratos.TIME] = self.analysis_start_time
        self.GetAnalysisModelPart().ProcessInfo[Kratos.STEP] = 0
        self.GetAnalysisModelPart().ProcessInfo[Kratos.IS_RESTARTED] = True
        self.analysis.Run()
            
    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return self.analysis._GetSolver().GetComputingModelPart()
    
    def GetAnalysisTimeSteppingData(self) -> 'tuple[float, float, Kratos.Parameters]':
        return (self.analysis_start_time, self.analysis_end_time, self.analysis_time_step_settings)

# test if execution policy works
if __name__ == "__main__":
    import json
    from matplotlib import pyplot as plt
    import KratosMultiphysics.StructuralMechanicsApplication as KratosSMA

    def CreateModelPart(model: Kratos.Model, model_part_name: str) -> None:
        model_part = model.CreateModelPart(model_part_name)

        model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(Kratos.VOLUME_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        model_part.AddNodalSolutionStepVariable(Kratos.POSITIVE_FACE_PRESSURE) 
        model_part.AddNodalSolutionStepVariable(Kratos.NEGATIVE_FACE_PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosSMA.POINT_LOAD)
        model_part.AddNodalSolutionStepVariable(KratosSMA.LINE_LOAD)
        model_part.AddNodalSolutionStepVariable(KratosSMA.SURFACE_LOAD)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)

        node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)

        prop = model_part.CreateNewProperties(1)
        element = model_part.CreateNewElement("TrussLinearElement3D2N", 1, [1, 2], prop)

        sub_model_part_A = model_part.CreateSubModelPart("support_A")
        sub_model_part_B = model_part.CreateSubModelPart("support_B")
        sub_model_part_Truss = model_part.CreateSubModelPart("truss")
        sub_model_part_A.AddNode(node1)
        sub_model_part_B.AddNode(node2)
        sub_model_part_Truss.AddNode(node1)
        sub_model_part_Truss.AddNode(node2)
        sub_model_part_Truss.AddElement(element) 

        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

    params = Kratos.Parameters("""{
            "model_part_names"        : ["Structure"],
            "analysis_module"         : "KratosMultiphysics.StructuralMechanicsApplication",
            "analysis_type"           : "StructuralMechanicsAnalysis",
            "analysis_settings"       : {
                "@include_json": "aux_files/project_parameters.json"
            }
    }""")

    model = Kratos.Model()

    # Construct and initialize execution policy
    CreateModelPart(model, "Structure")
    exec_policy = KratosTransientAnalysisExecutionPolicy("test", model, params)
    exec_policy.Initialize() # here the material file is read

    # First execution
    exec_policy.Execute()
    # exec_policy.Finalize()

    # Store result of first execution
    with open("aux_files/json_output_results.json", "r") as file:
        results = json.load(file)
    time_data = results["TIME"]
    disp_data = results["NODE_2"]["DISPLACEMENT_X"]

    # Second execution
    exec_policy.Execute()
    # exec_policy.Finalize()

    # Store result of second execution
    with open("aux_files/json_output_results.json", "r") as file:
        results = json.load(file)
    time_data2 = results["TIME"]
    disp_data2 = results["NODE_2"]["DISPLACEMENT_X"]

    # Change system parameters
    ##
    # Change system parameters using element expression
    # element_expression = Kratos.Expression.ElementExpression(model["Structure"])
    # Kratos.Expression.VariableExpressionIO.Read(element_expression, Kratos.YOUNG_MODULUS)
    # young_modulus_array = element_expression.Evaluate()
    # print(young_modulus_array)
    # young_modulus_array[0] *= 0.5  
    # Kratos.Expression.CArrayExpressionIO.Read(element_expression, young_modulus_array)
    # Kratos.Expression.VariableExpressionIO.Write(element_expression, Kratos.YOUNG_MODULUS)
    ##
    for index, element in enumerate(model["Structure"].Elements):
        element.Properties[Kratos.YOUNG_MODULUS] *= 0.5

    # Third execution
    exec_policy.Execute()
    # exec_policy.Finalize()

    # Store result of third execution
    with open("aux_files/json_output_results.json", "r") as file:
        results = json.load(file)
    time_data3 = results["TIME"]
    disp_data3 = results["NODE_2"]["DISPLACEMENT_X"]

    plt.plot(time_data, disp_data)
    plt.plot(time_data2, disp_data2, '.')
    plt.plot(time_data3, disp_data3)
    plt.show() 




        


