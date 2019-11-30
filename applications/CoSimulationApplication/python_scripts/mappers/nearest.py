from scipy.spatial import cKDTree
import numpy as np

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure

import time
from contextlib import contextmanager
@contextmanager
def timer(name=None, t=0, n=0, ms=False):
    startTime = time.time()
    yield
    elapsedTime = time.time() - startTime
    if ms:
        s = '\n' * n + '\t' * t + f'{elapsedTime * 1000:.2f}ms'
        s.replace(',', ' ')
    else:
        s = '\n' * n + '\t' * t + f'{elapsedTime:.1f}s'
    if name is not None:
        s += f' - {name}'
    s += '\n' * n
    print(s)

def Create(parameters):
    return MapperNearest(parameters)


# Class MapperNearest: 3D nearest-neighbour interpolation.
class MapperNearest(object):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters['settings']

        self.balanced_tree = self.settings['balanced_tree'].GetBool()

    def Initialize(self, model_part_from, model_part_to):
        coords_from = np.zeros((model_part_from.NumberOfNodes(), 3))
        for i, node in enumerate(model_part_from.Nodes):
            coords_from[i, :] = [node.X, node.Y, node.Z]

        coords_to = np.zeros((model_part_to.NumberOfNodes(), 3))
        for i, node in enumerate(model_part_to.Nodes):
            coords_to[i, :] = [node.X, node.Y, node.Z]

        # build and query tree
        if self.balanced_tree:  # time-intensive
            tree = cKDTree(coords_from)
        else:  # less stable
            tree = cKDTree(coords_from, balanced_tree=False)
        _, self.nearest = tree.query(coords_to)

    def Finalize(self):
        pass

    def __call__(self, args_from, args_to):
        model_part_from, var_from = args_from
        model_part_to, var_to = args_to

        # check if both Variables have same Type
        if var_from.Type() != var_to.Type():
            raise TypeError('Variables to be mapped have different Type.')

        # scalar interpolation
        if var_from.Type() == 'Double':
            hist_var_from = np.zeros(model_part_from.NumberOfNodes())
            for i, node in enumerate(model_part_from.Nodes):
                hist_var_from[i] = node.GetSolutionStepValue(var_from)

            hist_var_to = hist_var_from[self.nearest]  # nearest-neighbour interpolation
            for i, node in enumerate(model_part_to.Nodes):
                node.SetSolutionStepValue(var_to, 0, hist_var_to[i])

        # vector interpolation
        elif var_from.Type() == 'Array':
            hist_var_from = np.zeros((model_part_from.NumberOfNodes(), 3))
            for i, node in enumerate(model_part_from.Nodes):
                hist_var_from[i] = node.GetSolutionStepValue(var_from)

            hist_var_to = hist_var_from[self.nearest, :]  # nearest-neighbour interpolation
            for i, node in enumerate(model_part_to.Nodes):
                node.SetSolutionStepValue(var_to, 0, hist_var_to[i, :].tolist())

        # other types of Variables
        else:
            raise NotImplementedError(f'Mapping not yet implemented for Variable of Type {var_from.Type()}.')
