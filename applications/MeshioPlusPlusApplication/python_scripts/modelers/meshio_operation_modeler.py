import KratosMultiphysics
import KratosMultiphysics.MeshioPlusPlusApplication as KratosMeshioPlusPlus


class MeshioOperationModeler(KratosMultiphysics.Modeler):
    """Modeler applying a meshio++ mesh operation to a model part.

    The operation is selected with the "operation" setting, and its own settings are
    passed through in "operation_settings". Query what this build supports with
    KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.GetSupportedOperations(), and the
    available settings with GetDefaultParameters().

    The result is written into "output_model_part_name". The multi-output operations
    ("split", "partition") instead create one model part per piece, named
    "<output_model_part_name>_<operation>_<index>"; their names are listed in the report.
    """

    def __init__(self, model, settings):
        super().__init__(model, settings)

        settings.ValidateAndAssignDefaults(self.__GetDefaultSettings())

        self.model = model
        self.settings = settings

        input_model_part_name = settings["input_model_part_name"].GetString()
        if not input_model_part_name:
            raise Exception("Missing 'input_model_part_name' in input settings.")

        output_model_part_name = settings["output_model_part_name"].GetString()
        if not output_model_part_name:
            raise Exception("Missing 'output_model_part_name' in input settings.")

        # Created here, before the solvers add the nodal variables
        self.output_model_part = model.CreateModelPart(output_model_part_name)

    def SetupGeometryModel(self):
        super().SetupGeometryModel()

        input_model_part = self.model[self.settings["input_model_part_name"].GetString()]

        operation_settings = self.settings["operation_settings"].Clone()
        operation_settings.AddString("operation", self.settings["operation"].GetString())

        report = KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Execute(
            input_model_part, operation_settings, self.output_model_part)

        if self.settings["echo_level"].GetInt() > 0:
            KratosMultiphysics.Logger.PrintInfo(
                "MeshioOperationModeler",
                "Operation '{}': {}".format(self.settings["operation"].GetString(),
                                            report.PrettyPrintJsonString()))

    def PrepareGeometryModel(self):
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        super().SetupModelPart()

    @classmethod
    def __GetDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters("""{
            "echo_level"              : 0,
            "operation"               : "clean",
            "input_model_part_name"   : "",
            "output_model_part_name"  : "",
            "operation_settings"      : {}
        }""")
        return default_settings


def Factory(model, settings):
    return MeshioOperationModeler(model, settings)
