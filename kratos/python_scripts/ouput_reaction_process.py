import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IgaApplication


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return OuputReactionProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class OuputReactionProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "please_specify_model_part_name"
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]

    def ExecuteFinalizeSolutionStep(self):
        x_reaction = []
        y_reaction = []
        z_reaction = []

        for node in self.model_part.Nodes:
            reactions = node.GetSolutionStepValue(KratosMultiphysics.REACTION, 0)

            print(reactions)

            x_reaction.append(reactions[0])
            y_reaction.append(reactions[1])
            z_reaction.append(reactions[2])