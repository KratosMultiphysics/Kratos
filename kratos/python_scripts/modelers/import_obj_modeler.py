import KratosMultiphysics

def Factory(model, settings):
    """
    Factory function to create an instance of ImportOBJModeler.

    Parameters:
    model (KratosMultiphysics.Model): The model to which the OBJ model will be imported.
    settings (KratosMultiphysics.Parameters): The settings for the import process.

    Returns:
    ImportOBJModeler: An instance of ImportOBJModeler.
    """
    return ImportOBJModeler(model, settings)

class ImportOBJModeler(KratosMultiphysics.Modeler):
    """
    A class to import OBJ models into a KratosMultiphysics model part.

    Attributes:
    model_part (KratosMultiphysics.ModelPart): The model part where the OBJ model will be imported.
    settings (KratosMultiphysics.Parameters): The settings for the import process.

    Methods:
    SetupGeometryModel(): Sets up the geometry model.
    PrepareGeometryModel(): Prepares the geometry model.
    SetupModelPart(): Sets up the model part and imports the OBJ model.
    """

    def __init__(self, model, settings):
        """
        Initializes the ImportOBJModeler with the given model and settings.

        Parameters:
        model (KratosMultiphysics.Model): The model to which the OBJ model will be imported.
        settings (KratosMultiphysics.Parameters): The settings for the import process.
        """
        super().__init__(model, settings)

        # Cannot validate as settings may differ among input types
        settings.RecursivelyAddMissingParameters(self.__GetDefaultSettings())

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
        """
        Sets up the geometry model.
        """
        super().SetupGeometryModel()

    def PrepareGeometryModel(self):
        """
        Prepares the geometry model.
        """
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        """
        Sets up the model part and imports the OBJ model.

        Raises:
        Exception: If 'input_filename' is missing in the settings.
        """
        super().SetupModelPart()

        # Import the model part data
        # Note that at this point solvers must have already added the variables to the nodal variable data
        input_filename = self.settings["input_filename"].GetString()
        if not input_filename:
            err_msg = "Missing 'input_filename' in input settings."
            raise Exception(err_msg)

        # Ensure that NORMAL variable is added to the model part if required as historical
        obj_io_settings = self.settings["obj_io_settings"]
        if obj_io_settings["normal_as_historical"].GetBool() and not self.model_part.HasNodalSolutionStepVariable(KratosMultiphysics.NORMAL):
            self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        # Create ObjIO instance and read the model part
        obj_io = KratosMultiphysics.ObjIO(input_filename, obj_io_settings)
        obj_io.ReadModelPart(self.model_part)

    @classmethod
    def __GetDefaultSettings(cls):
        """
        Returns the default settings for the ImportOBJModeler.

        Returns:
        KratosMultiphysics.Parameters: The default settings.
        """
        default_settings = KratosMultiphysics.Parameters('''{
            "input_filename"  : "",
            "model_part_name" : "",
            "obj_io_settings"  : {
                "open_mode"                       : "read",
                "entity_type"                     : "element",
                "decompose_quads_into_triangles"  : false,
                "normal_as_historical"            : false
            }
        }''')
        return default_settings