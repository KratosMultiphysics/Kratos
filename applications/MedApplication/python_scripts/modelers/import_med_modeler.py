import KratosMultiphysics as KM
import KratosMultiphysics.MedApplication as KratosMed


class ImportMedModeler(KM.Modeler):
    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Cannot validate as settings may differ among input types
        settings.AddMissingParameters(self.__GetDefaultSettings())

        self.input_filename: str = settings["input_filename"].GetString()

        # Create the import destination model part
        # It is mandatory to do this when the modeler is instantiated to have the model part created before the solvers add the variables
        model_part_name: str = settings["model_part_name"].GetString()
        if not model_part_name:
            raise Exception(
                "Missing 'model_part_name' in input settings. This is where the imported model part is to be stored."
            )

        self.model_part: KM.ModelPart = model.CreateModelPart(model_part_name)

    def SetupGeometryModel(self):
        super().SetupGeometryModel()

        KratosMed.MedModelPartIO(self.input_filename, KM.IO.READ).ReadModelPart(self.model_part)

    def PrepareGeometryModel(self):
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        super().SetupModelPart()

    @classmethod
    def __GetDefaultSettings(cls):
        default_settings = KM.Parameters(
            """{
            "echo_level" : 0,
            "input_filename" : "",
            "model_part_name" : ""
        }"""
        )
        return default_settings


def Factory(model, settings):
    return ImportMedModeler(model, settings)
