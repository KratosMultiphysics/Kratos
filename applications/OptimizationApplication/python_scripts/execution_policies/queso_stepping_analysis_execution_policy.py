try:
    from QuESo_PythonApplication.PyQuESo import PyQuESo
except ImportError:
    raise Exception("QUESO python library is not available")

from importlib import import_module
import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos
import KratosMultiphysics.OptimizationApplication as KOA

def run_modelers(current_model, modelers_list):
    from KratosMultiphysics.modeler_factory import KratosModelerFactory
    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"QuesoSteppingAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"QuesoSteppingAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return QuesoSteppingAnalysisExecutionPolicy(parameters["name"].GetString(), model, parameters["settings"])

class QuesoSteppingAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "embedded_model_part_name" : "",
            "nurbs_model_part_name" : "",
            "analysis_module"  : "KratosMultiphysics",
            "analysis_type"    : "",
            "queso_json_input_file"    : "",
            "analysis_settings": {}
        }""")
        self.model = model
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)

        self.embedded_model_part = self.model.GetModelPart(parameters["embedded_model_part_name"].GetString())
        self.embedded_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        self.nurbs_model_part_name = parameters["nurbs_model_part_name"].GetString()

        self.analysis_module = parameters["analysis_module"].GetString()
        self.analysis_type = parameters["analysis_type"].GetString()
        self.analysis_settings = parameters["analysis_settings"]
        queso_json_input_file = parameters["queso_json_input_file"].GetString()

        if self.analysis_module == "KratosMultiphysics":
            self.analysis_full_module, self.analysis_type = GetClassModuleFromKratos(self.analysis_type)
        else:
            self.analysis_full_module = f"{self.analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(self.analysis_type)}"

        self.pyqueso = PyQuESo(queso_json_input_file)

        self.CreateAnalysis()

    def GetAnalysisModelPart(self):
        return self.analysis._GetSolver().GetComputingModelPart()

    def CreateAnalysis(self):
        hist_var_list = []
        if self.model.HasModelPart(self.nurbs_model_part_name):
            hist_var_list = self.model.GetModelPart(self.nurbs_model_part_name).GetHistoricalVariablesNames()
            self.model.DeleteModelPart(self.nurbs_model_part_name)
        self.model_part = self.model.CreateModelPart(self.nurbs_model_part_name)
        for var_name in hist_var_list:
            self.model_part.AddNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(var_name))
        self.analysis: AnalysisStage = getattr(import_module(self.analysis_full_module), self.analysis_type)(self.model, self.analysis_settings.Clone())

    def Initialize(self):

        # Read nurbs_volume_model_part
        modeler_settings = Kratos.Parameters("""
            [{
                "modeler_name": "NurbsGeometryModeler",
                "Parameters": {
                    "model_part_name" : "NurbsMesh1",
                    "geometry_name"   : "NurbsVolume"
                }
            }]
            """)
        queso_params = self.pyqueso.parameters
        tmp_parameters = modeler_settings[0]["Parameters"]
        tmp_parameters.AddEmptyValue("lower_point_xyz")
        tmp_parameters["lower_point_xyz"].SetVector(queso_params.LowerBoundXYZ())
        tmp_parameters.AddEmptyValue("upper_point_xyz")
        tmp_parameters["upper_point_xyz"].SetVector(queso_params.UpperBoundXYZ())
        tmp_parameters.AddEmptyValue("lower_point_uvw")
        tmp_parameters["lower_point_uvw"].SetVector(queso_params.LowerBoundUVW())
        tmp_parameters.AddEmptyValue("upper_point_uvw")
        tmp_parameters["upper_point_uvw"].SetVector(queso_params.UpperBoundUVW())
        tmp_parameters.AddEmptyValue("polynomial_order")
        tmp_parameters["polynomial_order"].SetVector(queso_params.Order())
        tmp_parameters.AddEmptyValue("number_of_knot_spans")
        tmp_parameters["number_of_knot_spans"].SetVector(queso_params.NumberOfElements())
        tmp_parameters["model_part_name"].SetString(self.nurbs_model_part_name)

        run_modelers(self.model, modeler_settings)
        self.pyqueso.Run(self.embedded_model_part)
        self.pyqueso.UpdateKratosNurbsVolumeModelPart(self.model_part)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def Execute(self):

        self.CreateAnalysis()

        self.Initialize()

        self.analysis.Run()



