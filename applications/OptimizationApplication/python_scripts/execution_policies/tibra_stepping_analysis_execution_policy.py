try:
    from TIBRA_PythonApplication.PyTIBRA import PyTIBRA
except ImportError:
    raise Exception("TIBRA python library is not available")

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
        raise RuntimeError(f"TibraSteppingAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"TibraSteppingAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return TibraSteppingAnalysisExecutionPolicy(parameters["name"].GetString(), model, parameters["settings"])

class TibraSteppingAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "embedded_model_part_name" : "",
            "nurbs_model_part_name" : "",
            "analysis_module"  : "KratosMultiphysics",
            "analysis_type"    : "",
            "analysis_settings": {}
        }""")
        self.model = model
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)

        self.embedded_model_part = self.model.GetModelPart(parameters["embedded_model_part_name"].GetString())
        self.nurbs_model_part_name = parameters["nurbs_model_part_name"].GetString()

        self.analysis_module = parameters["analysis_module"].GetString()
        self.analysis_type = parameters["analysis_type"].GetString()
        self.analysis_settings = parameters["analysis_settings"]

        if self.analysis_module == "KratosMultiphysics":
            self.analysis_full_module, self.analysis_type = GetClassModuleFromKratos(self.analysis_type)
        else:
            self.analysis_full_module = f"{self.analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(self.analysis_type)}"

        self.pytibra = PyTIBRA("TIBRAParameters.json")

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
                    "model_part_name" : "NurbsMesh",
                    "geometry_name"   : "NurbsVolume",
                    "lower_point": [-0.130, -0.110, -0.110],
                    "upper_point": [0.020, 0.190, 0.190],
                    "polynomial_order" : [2, 2, 2],
                    "number_of_knot_spans" : [50, 100, 100]
                }
            }]
            """)
        run_modelers(self.model, modeler_settings)
        self.pytibra.Clear()
        self.pytibra.Run(self.embedded_model_part)
        self.pytibra.UpdateKratosNurbsVolumeModelPart(self.model_part)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def Execute(self):

        self.CreateAnalysis()

        self.Initialize()

        self.analysis.Run()



