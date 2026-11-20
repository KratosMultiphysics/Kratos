# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ComputeBeamVectorsProcess(model, settings["Parameters"])

class ComputeBeamVectorsProcess(KratosMultiphysics.Process):
    """This process computes the tangential (T0) and normal (N0) vectors for isogeometric beam elements during initialization.
    """
    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""{
            "model_part_name" : ""
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.model = model
        self.model_part_beam = self.model[self.params["model_part_name"].GetString()]

        self.process = IGA.ComputeBeamVectorsProcess(self.model_part_beam)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()     