import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SM
import KratosMultiphysics.ConstitutiveLawsApplication as CLA

from matplotlib import pyplot as plt
import matplotlib.tri as tri
from matplotlib.pylab import *
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.animation as animation
import numpy as np

def Factory(settings, model):
    if(type(settings) != KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return OuputReactionProcess2(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class OuputReactionProcess2(KM.Process):
    def __init__(self, model, settings ):
        KM.Process.__init__(self)

        default_settings = KM.Parameters("""{
            "model_part_reactions" : "please_specify_model_part_name",
            "output_file_name" : "please_specify_output_file_name",
            "print_direction_index" : 0,
            "print_time_step": [0, 10]
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part_reactions = model[settings["model_part_reactions"].GetString()]
        self.model_part_root = self.model_part_reactions.GetRootModelPart()

        self.time_steps = settings["print_time_step"].GetVector()

        self.output_file_name_start = settings["output_file_name"].GetString()

        self.output_file_name = settings["output_file_name"].GetString() + str(self.model_part_root.NumberOfNodes())
        if not self.output_file_name.endswith(".txt"):
            self.output_file_name += ".txt"

        with open(self.output_file_name, 'w') as output_file:
            output_file.write("")

        self.print_direction_index = settings["print_direction_index"].GetInt()

        self.output_damage_file_name = "THE_DAMAGE_DATA_" + str(self.model_part_root.NumberOfNodes())
        if not self.output_damage_file_name.endswith(".txt"):
            self.output_damage_file_name += ".txt"

        with open(self.output_damage_file_name, 'w') as output_file:
            output_file.write("")

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

        self.line1, = self.axarr[0].plot(self.x_data, self.y_data_x_reactions, 'g')
        self.line2, = self.axarr[0].plot(self.x_data, self.y_data_y_reactions, 'b-')
        self.line3, = self.axarr[0].plot(self.x_data, self.y_data_z_reactions, 'b-')
        self.line4, = self.axarr[0].plot(self.x_data, self.y_data_length_reactions, 'r')



    def ExecuteFinalizeSolutionStep(self):
        time_step = self.model_part_root.ProcessInfo[KM.TIME]
        self.x_data.append(time_step)

        x_reaction = 0.0
        y_reaction = 0.0
        z_reaction = 0.0
        reaction = 0.0

        #for node in self.model_part_reactions.Nodes:
        #    reactions = node.GetSolutionStepValue(KM.REACTION, 0)

        #    x_reaction += reactions[0]
        #    y_reaction += reactions[1]
        #    z_reaction += reactions[2]

        #    reaction += sqrt(reactions[0]*reactions[0] + reactions[1]*reactions[1] + reactions[2]*reactions[2])

        for condition in self.model_part_reactions.Conditions:
            reactions = condition.CalculateOnIntegrationPoints(KM.REACTION, self.model_part_root.ProcessInfo)[0]

            x_reaction += reactions[0]
            y_reaction += reactions[1]
            z_reaction += reactions[2]

            reaction += sqrt(reactions[0]*reactions[0] + reactions[1]*reactions[1] + reactions[2]*reactions[2])

        with open(self.output_file_name, 'a') as output_file:
            #output_file.write(str(time_step)  + "  " + str(x_reaction) + "  " + str(y_reaction) + "  " + str(z_reaction) + "  " + str(reaction) + "\n")
            if self.print_direction_index == 0:
                output_file.write(str(time_step)  + "  " + str(x_reaction) + "\n")
            elif self.print_direction_index == 1:
                output_file.write(str(time_step)  + "  " + str(y_reaction) + "\n")
            elif self.print_direction_index == 2:
                output_file.write(str(time_step)  + "  " + str(z_reaction) + "\n")
            elif self.print_direction_index == 3:
                output_file.write(str(time_step)  + "  " + str(reaction) + "\n")

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


        x_location_elements = []
        y_location_elements = []
        z_location_elements = []

        x_reaction_elements = []
        y_reaction_elements = []
        z_reaction_elements = []

        solution_1 = []
        solution_2 = []

        for element in self.model_part_reactions.Elements:
            location = element.GetGeometry().Center()

            x_location_elements.append(location[0])
            y_location_elements.append(location[1])
            z_location_elements.append(location[2])

            N = element.GetGeometry().ShapeFunctionsValues()

            index = 0
            x_reaction = 0.0
            y_reaction = 0.0
            z_reaction = 0.0
            for node in element.GetGeometry():
                reactions = node.GetSolutionStepValue(KM.REACTION, 0)
                x_reaction = reactions[0] * N(0, index)
                y_reaction = reactions[1] * N(0, index)
                z_reaction = reactions[2] * N(0, index)

            x_reaction_elements.append(x_reaction)
            y_reaction_elements.append(y_reaction)
            z_reaction_elements.append(z_reaction)

            solution_1.append(0.0)
            solution_2.append(0.0)

        self.axarr[2].clear()
        self.axarr[2].quiver(x_location_elements, y_location_elements, x_reaction_elements, y_reaction_elements, scale=50)

        self.axarr[3].clear()
        self.axarr[3].scatter(x_location_elements, y_location_elements,s=10, c=solution_2, cmap='jet')

        plt.draw()
        #plt.savefig('reactions_' + str(time_step) + '.png')
        plt.pause(1e-17)

    def is_output_step(self, time_step):
        for t in self.time_steps:
            if (abs(time_step - t) < 10e-8):
                return True
        return False