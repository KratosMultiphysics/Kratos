import KratosMultiphysics as KM

class EmbeddedGeometryModeler(KM.Modeler):

    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Cannot validate as settings may differ ammong input types
        settings.AddMissingParameters(self.__GetDefaultSettings())

        # Declare required member variables
        self.model_part = None
        self.settings = settings

        # Get the model parts
        self.model_part = model.CreateModelPart(settings["model_part_name"].GetString())
        self.skin_model_part = model.CreateModelPart(settings["skin_model_part_name"].GetString())


    def SetupGeometryModel(self):
        super().SetupGeometryModel()

    def PrepareGeometryModel(self):
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        super().SetupModelPart()

        # Compute the level-set function
        domain_size = self.volume_model_part.ProcessInfo[KM.DOMAIN_SIZE]
        if domain_size == 2:
            if self.settings["discontinuous_distance"].GetBool():
                KM.CalculateDiscontinuousDistanceToSkinProcess2D(self.volume_model_part, self.embedded_model_part).Execute()
            else:
                KM.CalculateDistanceToSkinProcess2D(self.volume_model_part, self.embedded_model_part).Execute()
        elif domain_size == 3:
            if self.settings["discontinuous_distance"].GetBool():
                KM.CalculateDiscontinuousDistanceToSkinProcess3D(self.volume_model_part, self.embedded_model_part).Execute()
            else:
                KM.CalculateDistanceToSkinProcess3D(self.volume_model_part, self.embedded_model_part).Execute()
        else:
            raise Exception(f"Unknown DOMAIN_SIZE = {domain_size}")

        KM.Logger.PrintInfo("EmbeddedGeometryModeler", "Finished computing distance.")

    @classmethod
    def __GetDefaultSettings(cls):
        default_settings = KM.Parameters('''{
            "echo_level"               : 0,
            "volume_model_part_name"   : "",
            "embedded_model_part_name" : "",
            "discontinuous_distance"   : false
        }''')
        return default_settings

#TODO: Remove this once the modelers use the new model-parameters factory
def Factory(model, settings):
    return EmbeddedGeometryModeler(model, settings)
