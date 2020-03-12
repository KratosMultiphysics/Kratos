import numpy as np
import copy

from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperAxisymmetric3DTo2D(parameters)


# Class MapperAxisymmetric3DTo2D: map 3D to 2D axisymmetric.
class MapperAxisymmetric3DTo2D(CoSimulationComponent):
    def __init__(self, parameters):
        """
        - should there be both forward and backward initializations? NO
        - should there be a check to see whether input geometry is
            axisymmetrical wrt given directions?
        - take swirl into account? that would really be a
            3D axisymmetrical simulation (only geometry is 2D)...
        - TODO: documentation
        - TODO: create tests
        """
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

        if not forward:
            self.n_to = model_part_in.NumberOfNodes()


            self.n_from = self.n_t * self.n_to

            model = cs_data_structure.Model()
            model_part_out = model.CreateModelPart('no_name')
            model_part_out._ModelPart__hist_variables = model_part_in._ModelPart__hist_variables

            self.nearest = np.zeros((self.n_to, self.n_t)).astype(int)
            self.coeffs = np.ones((self.n_to, self.n_t)) / self.n_t
            self.theta = np.zeros((self.n_to, self.n_t))

            i_from = 0
            i_to = 0
            for node in model_part_in.Nodes:
                coords = np.array([node.X0, node.Y0, node.Z0])
                r = coords[self.dir_r]
                for i in range(self.n_t):
                    self.nearest[i_to, i] = i_from
                    self.theta[i_to, i] = i * 2 * np.pi / self.n_t

                    coords[self.dir_r] = r * np.cos(self.theta[i_to, i])
                    coords[self.dir_3d] = r * np.sin(self.theta[i_to, i])

                    model_part_out.CreateNewNode(i_from, *tuple(coords))
                    i_from += 1
                i_to += 1

            return model_part_out

        else:
            raise NotImplementedError('Forward Initialization not implemented for MapperAxisymmetric3DTo2D.')

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
                hist_var_to = np.dot(self.coeffs[i_to], hist_var_from[self.nearest[i_to]])
                node_to.SetSolutionStepValue(var_to, 0, hist_var_to)

        # vector interpolation
        elif var_from.Type() == 'Array':
            hist_var_from = np.zeros((self.n_from, 3))
            for i_from, node_from in enumerate(model_part_from.Nodes):
                hist_var_from[i_from] = node_from.GetSolutionStepValue(var_from)

            # *** old method: works
            # for i_to, node_to in enumerate(model_part_to.Nodes):
            #     v_a = 0.
            #     v_r = 0.
            #     for j in range(self.n_t):
            #         tmp = hist_var_from[self.nearest[i_to, j]]
            #         v_a += tmp[self.dir_a]
            #         theta = self.theta[i_to, j]
            #         v_r += tmp[self.dir_r] * np.cos(-theta) + tmp[self.dir_3d] * np.cos(np.pi / 2 - theta)
            #
            #     v_a /= self.n_t
            #     v_r /= self.n_t
            #     # v_r = 0.
            #
            #     hist_var_to = [0., 0., 0.]
            #     hist_var_to[self.dir_a] = v_a
            #     hist_var_to[self.dir_r] = v_r
            #     hist_var_to[self.dir_3d] = 0.
            #
            #     # print(hist_var_to)
            #
            #     node_to.SetSolutionStepValue(var_to, 0, hist_var_to)

            # *** new method: works too!
            for i_to, node_to in enumerate(model_part_to.Nodes):
                hist_var_to = [0., 0., 0.]
                tmp = hist_var_from[self.nearest[i_to]]
                hist_var_to[self.dir_a] = np.dot(self.coeffs[i_to], tmp[:, self.dir_a])
                hist_var_to[self.dir_r] = np.dot(self.coeffs[i_to],
                                                 tmp[:, self.dir_r] * np.cos(-self.theta[i_to]) +
                                                 tmp[:, self.dir_3d] * np.cos(np.pi / 2 - self.theta[i_to]))
                node_to.SetSolutionStepValue(var_to, 0, hist_var_to)

        # other types of Variables
        else:
            raise NotImplementedError(f'Mapping not yet implemented for Variable of Type {var_from.Type()}.')
