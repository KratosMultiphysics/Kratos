import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IgaApplication

from matplotlib import pyplot as plt
from matplotlib.pylab import *
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.animation as animation
import numpy

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return OuputReactionProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class OuputReactionProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "please_specify_model_part_name",
            "model_part_reactions" : "please_specify_model_part_name"
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.model_part_reactions = model[settings["model_part_reactions"].GetString()]

        self.f, self.axarr = plt.subplots(1, 4, figsize=(12, 4))

        self.x_data = []
        self.x_data.append(0.0)
        self.y_data_x_reactions = []
        self.y_data_x_reactions.append(0.0)
        self.y_data_y_reactions = []
        self.y_data_y_reactions.append(0.0)
        self.y_data_z_reactions = []
        self.y_data_z_reactions.append(0.0)
        self.y_data_length_reactions = []
        self.y_data_length_reactions.append(0.0)

        self.line1, = self.axarr[0].plot(self.x_data, self.y_data_x_reactions, 'b-')
        self.line2, = self.axarr[0].plot(self.x_data, self.y_data_y_reactions, 'b-')
        self.line3, = self.axarr[0].plot(self.x_data, self.y_data_z_reactions, 'b-')
        self.line4, = self.axarr[0].plot(self.x_data, self.y_data_length_reactions, 'b-')



    def ExecuteFinalizeSolutionStep(self):
        time_step = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.x_data.append(time_step)

        x_reaction = 0.0
        y_reaction = 0.0
        z_reaction = 0.0
        reaction = 0.0

        for node in self.model_part_reactions.Nodes:
            reactions = node.GetSolutionStepValue(KratosMultiphysics.REACTION, 0)

            x_reaction += reactions[0]
            y_reaction += reactions[1]
            z_reaction += reactions[2]

            reaction += sqrt(reactions[0]*reactions[0] + reactions[1]*reactions[1] + reactions[2]*reactions[2])

        self.y_data_x_reactions.append(x_reaction)
        self.y_data_y_reactions.append(y_reaction)
        self.y_data_z_reactions.append(z_reaction)
        self.y_data_length_reactions.append(reaction)

        self.line1.set_xdata(self.x_data)
        self.line1.set_ydata(self.y_data_x_reactions)
        self.line2.set_xdata(self.x_data)
        self.line2.set_ydata(self.y_data_y_reactions)
        self.line3.set_xdata(self.x_data)
        self.line3.set_ydata(self.y_data_z_reactions)
        self.line4.set_xdata(self.x_data)
        self.line4.set_ydata(self.y_data_length_reactions)

        self.axarr[0].relim()
        self.axarr[0].autoscale_view()

        x_location_nodes = []
        y_location_nodes = []
        z_location_nodes = []

        x_reaction_nodes = []
        y_reaction_nodes = []
        z_reaction_nodes = []

        for node in self.model_part.Nodes:
            reactions = node.GetSolutionStepValue(KratosMultiphysics.REACTION, 0)

            x_location_nodes.append(node.X)
            y_location_nodes.append(node.Y)
            z_location_nodes.append(node.Z)

            x_reaction_nodes.append(reactions[0])
            y_reaction_nodes.append(reactions[1])
            z_reaction_nodes.append(reactions[2])

        self.axarr[1].clear()
        self.axarr[1].quiver(x_location_nodes, y_location_nodes, x_reaction_nodes, y_reaction_nodes, scale=50)


        x_location_elements = []
        y_location_elements = []
        z_location_elements = []

        solution_1 = []
        solution_2 = []

        for element in self.model_part.Elements:
            location = element.Calculate(IgaApplication.COORDINATES, self.model_part.ProcessInfo)

            x_location_elements.append(location[0])
            y_location_elements.append(location[1])
            z_location_elements.append(location[2])

            solution_1.append(0.0)
            solution_2.append(0.0)

        self.axarr[2].clear()
        self.axarr[3].clear()

        self.axarr[2].scatter(x_location_elements, y_location_elements,s=10, c=solution_1, cmap='jet')
        self.axarr[3].scatter(x_location_elements, y_location_elements,s=10, c=solution_2, cmap='jet')

        plt.draw()
        plt.savefig('myfig' + str(time_step) + '.png')
        plt.pause(1e-17)
