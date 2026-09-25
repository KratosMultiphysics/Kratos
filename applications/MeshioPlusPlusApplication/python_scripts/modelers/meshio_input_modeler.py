import KratosMultiphysics
import KratosMultiphysics.MeshioPlusPlusApplication as KratosMeshioPlusPlus


class MeshioInputModeler(KratosMultiphysics.Modeler):
    """Modeler importing a model part from any meshio++-supported mesh format.

    The format is taken from the "input_format" setting, or resolved from the
    extension of "input_filename" when it is "auto". Query the formats
    available in this build with
    KratosMeshioPlusPlus.MeshioPlusPlusIO.GetSupportedReadFormats().

    "time_step" (default 0, negative counts from the end) selects one step of a
    multi-step file for the formats meshio++ reads selectively; "lenient"
    (default false) downgrades an unrepresentable mdpa/med construct to a
    warning instead of an error; "openfoam_region" (default "") selects one
    region of a multi-region OpenFOAM case. A partitioned file (vtkhdf, pvtu,
    pvtp, pvd, vtm) is read whole unless "select_piece" (default false) picks
    piece "piece" (default 0, negative counts from the end); "ghosts" ("keep"
    or "drop") removes the halo cells of a pvtu/pvtp/pvd. "mfem_grid_functions" and
    "patran_result_files" ({"name", "path"} lists) read fields stored next to an
    MFEM/Patran mesh; "z88_results" : false skips a Z88 deck's results;
    "read_field_data" carries the file's data onto the matching registered
    variables. See
    MeshioPlusPlusIO.ReadModelPart.
    """

    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Cannot validate as settings may differ among input types
        settings.AddMissingParameters(self.__GetDefaultSettings())

        # Declare required member variables
        self.model_part = None
        self.settings = settings

        # Create the import destination model part
        # It is mandatory to do this when the modeler is instantiated to have the model part created before the solvers add the variables
        model_part_name = self.settings["model_part_name"].GetString()
        if not model_part_name:
            err_msg = "Missing 'model_part_name' in input settings. This is where the imported model part is to be stored."
            raise Exception(err_msg)
        else:
            self.model_part = model.CreateModelPart(model_part_name)

    def SetupGeometryModel(self):
        super().SetupGeometryModel()

        # Import the model part data
        # Note that at this point solvers must have already added the variables to the nodal variable data
        io_settings = KratosMultiphysics.Parameters("""{}""")
        io_settings.AddString("format", self.settings["input_format"].GetString())
        io_settings.AddInt("time_step", self.settings["time_step"].GetInt())
        io_settings.AddBool("lenient", self.settings["lenient"].GetBool())
        io_settings.AddString("openfoam_region", self.settings["openfoam_region"].GetString())
        io_settings.AddBool("select_piece", self.settings["select_piece"].GetBool())
        io_settings.AddInt("piece", self.settings["piece"].GetInt())
        io_settings.AddString("ghosts", self.settings["ghosts"].GetString())
        io_settings.AddValue("mfem_grid_functions", self.settings["mfem_grid_functions"])
        io_settings.AddValue("patran_result_files", self.settings["patran_result_files"])
        io_settings.AddBool("z88_results", self.settings["z88_results"].GetBool())
        io_settings.AddBool("read_field_data", self.settings["read_field_data"].GetBool())
        meshio_io = KratosMeshioPlusPlus.MeshioPlusPlusIO(
            self.settings["input_filename"].GetString(),
            io_settings)
        meshio_io.ReadModelPart(self.model_part)

    def PrepareGeometryModel(self):
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        super().SetupModelPart()

    @classmethod
    def __GetDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters('''{
            "echo_level"       : 0,
            "input_filename"   : "",
            "input_format"     : "auto",
            "model_part_name"  : "",
            "time_step"        : 0,
            "lenient"          : false,
            "openfoam_region"  : "",
            "select_piece"     : false,
            "piece"            : 0,
            "ghosts"           : "keep",
            "mfem_grid_functions" : [],
            "patran_result_files" : [],
            "z88_results"      : true,
            "read_field_data"  : false
        }''')
        return default_settings


def Factory(model, settings):
    return MeshioInputModeler(model, settings)
