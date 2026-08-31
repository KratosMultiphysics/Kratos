import KratosMultiphysics
import KratosMultiphysics.MeshioPlusPlusApplication as KratosMeshioPlusPlus


class MeshioOperationModeler(KratosMultiphysics.Modeler):
    """Modeler applying a meshio++ mesh operation to a model part.

    The operation is selected with the "operation" setting, and its own settings are
    passed through in "operation_settings" - including "clean", "transform", "convert_cells",
    "refine", "decimate", "smooth", "reorder", "extract_surface", "extract_skin", "crop_bbox",
    "crop_halfspace", "crop_predicate", "slice", "isosurface", "attach_quality", "gradient",
    "voxelize", "compute_sdf", "split", "partition", "stats", "quality", and the data
    operations "data_calc", "data_condition", "data_manage", "data_info",
    "point_data_to_cell_data", "cell_data_to_point_data". Query what this build supports with
    KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.GetSupportedOperations(), and the
    available settings with GetDefaultParameters(). "interpolate" is not reachable here - it
    needs two source model parts; see MeshioInterpolateModeler. Nor are "Grid" (no source at
    all), "DistanceToSurface" (two sources) and "CheckSurfaceWatertight" (report-only); call
    those on MeshioPlusPlusMeshOperations directly.

    Field data (nodal/elemental/conditional variables, flags and ids) can be carried through
    the operation with the same "nodal_solution_step_data_variables" /
    "nodal_data_value_variables" / "nodal_flags" / "element_data_value_variables" /
    "element_flags" / "condition_data_value_variables" / "condition_flags" /
    "gauss_point_variables_in_elements" / "write_ids" / "write_mdpa_ids" settings
    MeshioPlusPlusIO uses, empty by default. A resulting array is written back onto the output
    model part only when it is Float64 or Int64 and its name matches a registered Variable with
    the right component count - point "data_calc"'s "output" setting (or "gradient"'s,
    "compute_sdf"'s, ...) at an existing variable name to get the result back.

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
