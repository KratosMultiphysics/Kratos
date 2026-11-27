import KratosMultiphysics
import KratosMultiphysics.MeshioPlusPlusApplication as KratosMeshioPlusPlus


class MeshioInterpolateModeler(KratosMultiphysics.Modeler):
    """Modeler sampling one model part's field data onto another's geometry.

    Unlike MeshioOperationModeler this needs two independent sources - "source_model_part_name"
    (carrying the field data to sample) and "target_model_part_name" (contributing only its
    geometry) - so it is a separate modeler rather than another "operation" entry.

    Field data (nodal/elemental/conditional variables, flags and ids) is selected with the
    same settings MeshioPlusPlusIO uses ("nodal_solution_step_data_variables",
    "nodal_data_value_variables", "nodal_flags", "element_data_value_variables",
    "element_flags", "condition_data_value_variables", "condition_flags",
    "gauss_point_variables_in_elements", "write_ids"), applied to *both* the source and the
    target - the usual case is that only the source carries the named data, but a name present
    on both is exactly what "on_conflict" resolves. Interpolation itself is configured through
    "operation_settings" ("method": "nearest"/"barycentric"; "names": array names to
    interpolate, empty = every array the source's field data settings collect; "extrapolate";
    "default_value"; "on_conflict": "error"/"overwrite"/"suffix"). Query the available
    settings with KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.GetDefaultParameters().

    The result - the target's topology plus the interpolated arrays - is written into
    "output_model_part_name".
    """

    def __init__(self, model, settings):
        super().__init__(model, settings)

        settings.ValidateAndAssignDefaults(self.__GetDefaultSettings())

        self.model = model
        self.settings = settings

        source_model_part_name = settings["source_model_part_name"].GetString()
        if not source_model_part_name:
            raise Exception("Missing 'source_model_part_name' in input settings.")

        target_model_part_name = settings["target_model_part_name"].GetString()
        if not target_model_part_name:
            raise Exception("Missing 'target_model_part_name' in input settings.")

        output_model_part_name = settings["output_model_part_name"].GetString()
        if not output_model_part_name:
            raise Exception("Missing 'output_model_part_name' in input settings.")

        # Created here, before the solvers add the nodal variables
        self.output_model_part = model.CreateModelPart(output_model_part_name)

    def SetupGeometryModel(self):
        super().SetupGeometryModel()

        source_model_part = self.model[self.settings["source_model_part_name"].GetString()]
        target_model_part = self.model[self.settings["target_model_part_name"].GetString()]

        report = KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Interpolate(
            source_model_part, target_model_part, self.settings["operation_settings"],
            self.output_model_part)

        if self.settings["echo_level"].GetInt() > 0:
            KratosMultiphysics.Logger.PrintInfo(
                "MeshioInterpolateModeler", report.PrettyPrintJsonString())

    def PrepareGeometryModel(self):
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        super().SetupModelPart()

    @classmethod
    def __GetDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters("""{
            "echo_level"              : 0,
            "source_model_part_name"  : "",
            "target_model_part_name"  : "",
            "output_model_part_name"  : "",
            "operation_settings"      : {}
        }""")
        return default_settings


def Factory(model, settings):
    return MeshioInterpolateModeler(model, settings)
