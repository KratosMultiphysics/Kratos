import KratosMultiphysics

class ImportMDPAModeler(KratosMultiphysics.Modeler):

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
        input_type = "mdpa"
        KratosMultiphysics.SingleImportModelPart.Import(
            self.model_part,
            self.settings,
            input_type)

    def PrepareGeometryModel(self):
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        super().SetupModelPart()

    @classmethod
    def __GetDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters('''{
            "echo_level" : 0,
            "input_filename" : "",
            "model_part_name" : ""
        }''')
        return default_settings

def Factory(model, settings):
    return ImportMDPAModeler(model, settings)
