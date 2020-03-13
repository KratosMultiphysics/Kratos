import numpy as np
import copy

from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperAxisymmetric2DTo3D(parameters)


# Class MapperAxisymmetric2DTo3D: map 2D axisymmetric to 3D.
class MapperAxisymmetric2DTo3D(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters['settings']
        self.interpolator = False

        # get axial and radial directions
        dirs = ['X', 'Y', 'Z']
        self.dir_a = dirs.index(self.settings['direction_axial'].GetString())
        self.dir_r = dirs.index(self.settings['direction_radial'].GetString())
        self.dir_3d = (set([0, 1 ,2]) - set([self.dir_a, self.dir_r])).pop()

        # get number of nodes in tangential direction
        self.n_t = self.settings['n_tangential'].GetInt()  # int parameter
        if self.n_t < 6:  # hardcode minimum value?
            raise ValueError('minimum value for n_tangential is 6')

    def Initialize(self, model_part_in, forward):
        super().Initialize()

        if forward:
            self.n_from = model_part_in.NumberOfNodes()

            n_t = self.settings['n_tangential'].GetInt()  # int parameter
            if n_t < 6:  # hardcode minimum value?
                raise ValueError('minimum value for n_tangential is 6')
            self.n_to = n_t * self.n_from

            model = cs_data_structure.Model()
            model_part_out = model.CreateModelPart('no_name')
            model_part_out._ModelPart__hist_variables = model_part_in._ModelPart__hist_variables

            self.nearest = np.zeros((self.n_to, 1)).astype(int)
            self.theta = np.zeros((self.n_to, 1))

            i_from = 0
            i_to = 0
            for node in model_part_in.Nodes:
                coords = np.array([node.X0, node.Y0, node.Z0])
                r = coords[self.dir_r]
                for i in range(n_t):
                    self.nearest[i_to, 0] = i_from
                    self.theta[i_to] = i * 2 * np.pi / n_t

                    coords[self.dir_r] = r * np.cos(self.theta[i_to])
                    coords[self.dir_3d] = r * np.sin(self.theta[i_to])

                    model_part_out.CreateNewNode(i_to, *tuple(coords))
                    i_to += 1
                i_from += 1

            return model_part_out

        else:
            raise NotImplementedError('Backward Initialization not implemented for MapperAxisymmetric2DTo3D.')


    def __call__(self, args_from, args_to):
        model_part_from, var_from = args_from
        model_part_to, var_to = args_to

        # check if both Variables have same Type
        if var_from.Type() != var_to.Type():
            raise TypeError('Variables to be mapped have different Type.')

        # scalar interpolation
        if var_from.Type() == 'Double':
            hist_var_from = np.zeros(self.n_from)
            for i_from, node_from in enumerate(model_part_from.Nodes):
                hist_var_from[i_from] = node_from.GetSolutionStepValue(var_from)

            for i_to, node_to in enumerate(model_part_to.Nodes):
                node_to.SetSolutionStepValue(var_to, 0, hist_var_from[self.nearest[i_to, 0]])

        # vector interpolation
        elif var_from.Type() == 'Array':
            hist_var_from = np.zeros((self.n_from, 3))
            for i_from, node_from in enumerate(model_part_from.Nodes):
                hist_var_from[i_from] = node_from.GetSolutionStepValue(var_from)

            for i_to, node_to in enumerate(model_part_to.Nodes):
                hist_var_to = hist_var_from[self.nearest[i_to, 0]].copy()

                tmp = hist_var_to[self.dir_r]
                hist_var_to[self.dir_r] = tmp * np.cos(self.theta[i_to])
                hist_var_to[self.dir_3d] = tmp * np.sin(self.theta[i_to])

                node_to.SetSolutionStepValue(var_to, 0, hist_var_to)

        # other types of Variables
        else:
            raise NotImplementedError(f'Mapping not yet implemented for Variable of Type {var_from.Type()}.')
