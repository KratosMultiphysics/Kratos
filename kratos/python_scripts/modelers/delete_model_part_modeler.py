import KratosMultiphysics

class DeleteModelPartModeler(KratosMultiphysics.Modeler):

    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Cannot validate as settings may differ ammong input types
        settings.AddMissingParameters(self.__GetDefaultSettings())

        # Declare required member variables
        self.model = model
        self.settings = settings

    def SetupGeometryModel(self):
        super().SetupGeometryModel()

    def PrepareGeometryModel(self):
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        super().SetupModelPart()

        # Remove the specified model part
        model_part_name = self.settings["model_part_name"].GetString()
        self.model.DeleteModelPart(model_part_name)
        if self.settings["echo_level"].GetInt() > 0:
            KratosMultiphysics.Logger.PrintInfo("DeleteModelPartModeler", f"Model part '{model_part_name}' deleted.")

    def __GetDefaultSettings(self):
        default_settings = KratosMultiphysics.Parameters('''{
            "echo_level" : 0,
            "model_part_name" : ""
        }''')
        return default_settings

def Factory(model, settings):
    return DeleteModelPartModeler(model, settings)