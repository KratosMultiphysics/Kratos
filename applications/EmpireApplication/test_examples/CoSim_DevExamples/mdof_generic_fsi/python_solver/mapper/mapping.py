# Class for the fluid part
# import numerical tools

import numpy as np
from math import sin, cos, radians
from KratosMultiphysics import *
from scipy.spatial import distance


class Mapper:

    def __init__(self, model_part, structure):

        self.model_part = model_part.Nodes
        self.structure = structure
        self.nodes = self.sort_nodes()
        self.forces = None
        self.mapped_forces = None

    # def sort_nodes(self):

    #     num_levels = self.structure.properties.levels
    #     level_height = self.structure.properties.height / num_levels

    #     nodes = [None] * num_levels

    #     for i in range(0, num_levels):
    #         nodes[i] = []
    #         for node in self.model_part:
    #             if i == 0:
    #                 if (i * level_height < node.Z0) and (node.Z0 <= (i + 1) * level_height):
    #                     nodes[i].append(node)
    #             else:
    #                 if (i * level_height < node.Z0) and (node.Z0 <= (i + 1) * level_height):
    #                     nodes[i].append(node)

    #     return nodes

    def sort_nodes(self):

        num_levels = self.structure.properties.levels
        level_height = self.structure.properties.height / num_levels

        struct_levels = {}

        for i in range(0, num_levels):
            struct_levels["level_" + str(i)] = []
            for node in self.model_part:
                if i == 0:
                    if (i * level_height <= node.Z0) and (node.Z0 <= (i + 1) * level_height):
                        struct_levels["level_" + str(i)].append(node)
                else:
                    if (i * level_height < node.Z0) and (node.Z0 <= (i + 1) * level_height):
                        struct_levels["level_" + str(i)].append(node)

        return struct_levels

    # Extract forces from levels
    # def extract_forces(self):
    #     c = 0
    #     num_levels = self.structure.properties.levels
    #     F_X, F_Y, M_Z = np.zeros(num_levels), np.zeros(
    #         num_levels), np.zeros(num_levels)
    #     f_x = 0

    #     for level in range(len(self.nodes)):
    #         level_reaction = [0.0, 0.0, 0.0]
    #         level_moment = 0

    #         theta = radians(self.structure.results[5][level])
    #         print("THETA STRUCTURE", theta)
    # input()
    #         T = np.array([[cos(theta), sin(theta), 0.],
    #                       [sin(theta), cos(theta), 0.],
    #                       [0., 0., 1.]])

    #         for node in self.nodes[level]:
    #             node_position = [node.X, node.Y, node.Z]

    #             pos_vector = [a - b for a, b in zip(
    #                 node_position, self.structure.position[level])]

    #             nodal_result = node.GetSolutionStepValue(REACTION, 0)
    #             nodal_result_rot = np.dot(T, nodal_result)
    #             level_reaction = level_reaction - nodal_result_rot

    #             level_moment = level_moment + \
    #                 (-nodal_result_rot[0] * pos_vector[1] +
    #                  nodal_result_rot[1] * pos_vector[0])
    #             f_x += node.GetSolutionStepValue(REACTION, 0)[0]
    #             c+=1
    #         level_reaction[2] = level_moment

    #         F_X[level] = level_reaction[0]
    #         F_Y[level] = level_reaction[1]
    #         M_Z[level] = level_reaction[2]
    #     print("F_X: ", f_x)
    #     print("NR NODES", c)
    #     print("EXTRACTED FORCES:", F_X)
    #     input()
    #     self.forces = F_X, F_Y, M_Z

    def extract_forces(self):
        num_levels = self.structure.properties.levels
        F_X, F_Y, M_Z = np.zeros(num_levels), np.zeros(
            num_levels), np.zeros(num_levels)
        f_x = 0
        l = 0

        for level, nodes in sorted(self.nodes.items()):
            level_reaction = [0.0, 0.0, 0.0]
            level_moment = 0

            theta = radians(self.structure.results[5][l])

            T = np.array([[cos(theta), sin(theta), 0.],
                          [sin(theta), cos(theta), 0.],
                          [0., 0., 1.]])

            for node in nodes:
                node_position = [node.X, node.Y, node.Z]

                pos_vector = [a - b for a, b in zip(
                    node_position, self.structure.position[l])]

                nodal_result = node.GetSolutionStepValue(REACTION, 0)
                nodal_result_rot = np.dot(T, nodal_result)
                level_reaction = level_reaction - nodal_result_rot

                level_moment = level_moment + \
                    (-nodal_result_rot[0] * pos_vector[1] +
                     nodal_result_rot[1] * pos_vector[0])
                f_x += node.GetSolutionStepValue(REACTION, 0)[0]

            level_reaction[2] = level_moment

            F_X[l] = level_reaction[0]
            F_Y[l] = level_reaction[1]
            M_Z[l] = level_reaction[2]
            l += 1

        self.forces = F_X, F_Y, M_Z

    def map_forces_to_structure(self):

        def load_distribution(extracted_forces):
            level_number = self.structure.properties.levels
            level_length = self.structure.properties.height / level_number

            def nodal_force(force):
                F = [force / 2, force * level_length / 12,
                     force / 2, -force * level_length / 12]
                return F

            nodal_load = list(map(nodal_force, extracted_forces))

            load = np.zeros(2 * level_number + 2)

            for i in range(level_number):
                load_temp = np.zeros(2 * level_number + 2)
                load_temp[2 * i:2 * i + 4] = nodal_load[i]
                load += load_temp

            # remove the fixed degrees of freedom
            rdof = [1, 0]
            for dof in rdof:
                load = np.delete(load, dof, axis=0)

            return load

        def load_distribution_torsion(extracted_forces):
            level_number = self.structure.properties.levels
            level_length = self.structure.properties.height / level_number

            def nodal_force(force):
                F = [force / 2, force / 2]
                return F

            nodal_load = list(map(nodal_force, extracted_forces))

            load = np.zeros(level_number + 1)

            for i in range(level_number):
                load_temp = np.zeros(level_number + 1)
                load_temp[i:i + 2] = nodal_load[i]
                load += load_temp

            # remove the fixed degrees of freedom
            rdof = [0]
            for dof in rdof:
                load = np.delete(load, dof, axis=0)

            return load

        mapped_force_X = load_distribution(self.forces[0])
        mapped_force_Y = load_distribution(self.forces[1])
        mapped_force_R = load_distribution_torsion(self.forces[2])

        self.mapped_forces = mapped_force_X, mapped_force_Y, mapped_force_R

    def map_from_file_to_structure(self, from_file):

        def load_distribution(forces):
            level_number = self.structure.properties.levels
            level_length = self.structure.properties.height / level_number

            def nodal_force(force):
                F = [force / 2, force * level_length / 12,
                     force / 2, -force * level_length / 12]
                return F

            nodal_load = list(map(nodal_force, forces))

            load = np.zeros(2 * level_number + 2)

            for i in range(level_number):
                load_temp = np.zeros(2 * level_number + 2)
                load_temp[2 * i:2 * i + 4] = nodal_load[i]
                load += load_temp

            # remove the fixed degrees of freedom
            rdof = [1, 0]
            for dof in rdof:
                load = np.delete(load, dof, axis=0)

            return load

        mapped_force_X = load_distribution(from_file[0])
        mapped_force_Y = load_distribution(from_file[1])
        mapped_force_R = from_file[2]

        return [mapped_force_X, mapped_force_Y, mapped_force_R]

    def nodal_displacements(self, res, level, node):

        results = [None] * len(res)
        for i in range(len(res)):
            results[i] = np.zeros(len(res[i]) + 1)
            results[i][1:] = res[i]

        num_levels = self.structure.properties.levels
        level_length = self.structure.properties.height / num_levels

        height = (level + 1) * level_length
        node_rel_pos = height - node.Z

        xi = -(node_rel_pos / level_length) + 1
        alpha = xi * (
            results[5][level + 1] - results[5][level]) + results[5][level]
        beta = xi * (
            results[4][level + 1] - results[4][level]) + results[4][level]
        gamma = xi * (
            results[3][level + 1] - results[3][level]) + results[3][level]
        disp_x = xi * (
            results[0][level + 1] - results[0][level]) + results[0][level]
        disp_y = xi * (
            results[1][level + 1] - results[1][level]) + results[1][level]
        disp_z = xi * (
            results[2][level + 1] - results[2][level]) + results[2][level]

        return [alpha, beta, gamma, disp_x, disp_y, disp_z]

    def transformation_matrix(self, nodal_values):

        alpha, beta, gamma, dispX, dispY, dispZ = nodal_values

        # Transformation Matrix
        T = np.zeros((4, 4))
        T[0, 0] = cos(radians(alpha)) * cos(radians(beta))
        T[0, 1] = cos(radians(alpha)) * sin(radians(beta)) * sin(
            radians(gamma)) - sin(radians(alpha)) * cos(radians(gamma))
        T[0, 2] = cos(radians(alpha)) * sin(radians(beta)) * cos(
            radians(gamma)) + sin(radians(alpha)) * sin(radians(gamma))
        T[0, 3] = dispX
        T[1, 0] = sin(radians(alpha)) * cos(radians(beta))
        T[1, 1] = sin(radians(alpha)) * sin(radians(beta)) * sin(
            radians(gamma)) + cos(radians(alpha)) * cos(radians(gamma))
        T[1, 2] = sin(radians(alpha)) * sin(radians(beta)) * cos(
            radians(gamma)) - cos(radians(alpha)) * sin(radians(gamma))
        T[1, 3] = dispY
        T[2, 0] = -sin(radians(beta))
        T[2, 1] = cos(radians(beta)) * sin(radians(gamma))
        T[2, 2] = cos(radians(beta)) * cos(radians(gamma))
        T[2, 3] = dispZ
        T[3, 0] = 0
        T[3, 1] = 0
        T[3, 2] = 0
        T[3, 3] = 1
        return T

    def set_mesh_displacement(self):
        l = 0
        for level, nodes in sorted(self.nodes.items()):

            for node in nodes:
                nodal_values = self.nodal_displacements(
                    self.structure.results, l, node)

                # Transformation Matrix
                T = self.transformation_matrix(nodal_values)

                r_0 = [node.X - 0.0, node.Y - 0.0, node.Z - 0.0, 1]
                r = np.dot(T, r_0)

                # Elementary displacements
                dx = r[0] - r_0[0]
                dy = r[1] - r_0[1]
                dz = r[2] - r_0[2]

                # Set solution to mesh
                node.SetSolutionStepValue(MESH_DISPLACEMENT_X, dx)
                node.SetSolutionStepValue(MESH_DISPLACEMENT_Y, dy)
                node.SetSolutionStepValue(MESH_DISPLACEMENT_Z, dz)
            l += 1

                # if node.Id == 624:
                #     print("Nodal Values:", nodal_values)
                #     print("Transformation Matrix:", T)
                #     print("R_0:", r_0)
                #     print("R:", r)
                #     print("DISPLACED VALUE", dy)
                #     v1 = [node.X + dx, node.Y + dy, node.Z + dz]
                #     print(v1)
                # elif node.Id == 426:
                #     print("Nodal Values:", nodal_values)
                #     print("Transformation Matrix:", T)
                #     print("R_0:", r_0)
                #     print("R:", r)
                #     print("DISPLACED VALUE", dy)
                #     v2 = [node.X + dx, node.Y + dy, node.Z + dz]
        # print("VECTOR 1:", v1)
        # print("VECTOR 2:", v2)
        # print("Distance:", distance.euclidean(v1,v2))
        # print("Z displacement:", v2[2] - v1[2])
    def set_mesh_velocity_to_fluid(self):
        for node in self.model_part:
            # get mesh velocity at current step
            MeshVelocity = node.GetSolutionStepValue(MESH_VELOCITY, 0)
            # assign it to fluid velocity at current step
            node.SetSolutionStepValue(VELOCITY, 0, MeshVelocity)
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.Fix(VELOCITY_Z)
