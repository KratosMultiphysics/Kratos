import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SM

from matplotlib import pyplot as plt
import matplotlib.tri as tri
from matplotlib.pylab import *
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.animation as animation
import numpy as np

def Factory(settings, model):
    if(type(settings) != KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return OuputReactionProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class OuputReactionProcess(KM.Process):
    def __init__(self, model, settings ):
        KM.Process.__init__(self)

        default_settings = KM.Parameters("""{
            "model_part_name" : "please_specify_model_part_name",
            "output_file_name" : "reactions"
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]

        self.output_file_name = settings["output_file_name"].GetString()

        self.output_file_name_total = settings["output_file_name"].GetString() + "_total"
        if not self.output_file_name_total.endswith(".txt"):
            self.output_file_name_total += ".txt"

        with open(self.output_file_name_total, 'w') as output_file:
            output_file.write("")

        self.time_data = []

    def ExecuteFinalizeSolutionStep(self):
        time_step = self.model_part.ProcessInfo[KM.TIME]
        self.time_data.append(time_step)

        x_reaction = 0.0
        y_reaction = 0.0
        z_reaction = 0.0
        reaction = 0.0

        with open(self.output_file_name + "_" + str(time_step) + ".txt", 'w') as output_file:
            output_file.write("")

        for condition in self.model_part.Conditions:
            reactions = condition.CalculateOnIntegrationPoints(IGA.PENALTY_REACTION, self.model_part.ProcessInfo)[0]
            location = condition.GetGeometry().Center()

            with open(self.output_file_name + "_" + str(time_step) + ".txt", 'a') as output_file:
                output_file.write(str(location[0])  + "  " + str(location[1]) + "  " + str(location[2]) + "  " + str(reactions[0]) + "  " + str(reactions[1]) + "  " + str(reactions[2]) + "\n")

            x_reaction += reactions[0]
            y_reaction += reactions[1]
            z_reaction += reactions[2]

            reaction += sqrt(reactions[0]*reactions[0] + reactions[1]*reactions[1] + reactions[2]*reactions[2])

        with open(self.output_file_name_total, 'a') as output_file:
            output_file.write(str(time_step)  + "  " + str(x_reaction) + "  " + str(y_reaction) + "  " + str(z_reaction) + "  " + str(reaction) + "\n")