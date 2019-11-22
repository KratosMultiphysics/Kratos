from scipy.spatial import cKDTree
import numpy as np

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperNearest(parameters)


# Class MapperNearest: 3D nearest-neighbour interpolation.
class MapperNearest(object):
    def __init__(self, _unused):
        super().__init__()

    def Initialize(self, model_part_from, model_part_to):
        coords_from = np.zeros((model_part_from.NumberOfNodes(), 3))
        for i, node in enumerate(model_part_from.Nodes):
            coords_from[i, :] = [node.X, node.Y, node.Z]

        coords_to = np.zeros((model_part_to.NumberOfNodes(), 3))
        for i, node in enumerate(model_part_to.Nodes):
            coords_to[i, :] = [node.X, node.Y, node.Z]

        # find nearest neighbour
        tree = cKDTree(coords_from)
        _, self.nearest = tree.query(coords_to, n_jobs=-1)  # runs in parallel


        # *** add support for 3D variables
        # *** add check to see if variables have same dimensions

    def Finalize(self):
        pass

    def __call__(self, args_from, args_to):
        model_part_from, var_from = args_from
        model_part_to, var_to = args_to

        hist_var_from = np.zeros(model_part_from.NumberOfNodes())
        for i, node in enumerate(model_part_from.Nodes):
            hist_var_from[i] = node.GetSolutionStepValue(var_from)

        hist_var_to = hist_var_from[self.nearest]  # nearest-neighbour interpolation
        for i, node in enumerate(model_part_to.Nodes):
            node.SetSolutionStepValue(var_to, 0, hist_var_to[i])
