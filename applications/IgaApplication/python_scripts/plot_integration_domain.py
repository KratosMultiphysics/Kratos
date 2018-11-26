import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IgaApplication

from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import random

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PlotIntegrationDomain(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class PlotIntegrationDomain(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "please_specify_model_part_name",
            "store_plots"     : true
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        fig = pyplot.figure()
        ax = Axes3D(fig)

        x_coords = []
        y_coords = []
        z_coords = []

        for element in self.model_part.Elements:
            coords = element.Calculate(IgaApplication.COORDINATES, self.model_part.ProcessInfo)
            x_coords.append(coords[0])
            y_coords.append(coords[1])
            z_coords.append(coords[2])

        x_coords_nodes = []
        y_coords_nodes = []
        z_coords_nodes = []

        for node in self.model_part.GetRootModelPart().Nodes:
            x_coords_nodes.append(node.X)
            y_coords_nodes.append(node.Y)
            z_coords_nodes.append(node.Z)

        ax.scatter(x_coords, y_coords, z_coords)
        ax.scatter(x_coords_nodes, y_coords_nodes, z_coords_nodes)
        pyplot.show()
