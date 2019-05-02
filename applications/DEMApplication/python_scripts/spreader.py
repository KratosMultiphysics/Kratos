import math
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

def GetPolarRCoordinate2(node):
    return node.X ** 2 + node.Y ** 2

def GetPolarCoordinates(node):
    x = math.sqrt(node.X ** 2 + node.Y ** 2)
    y = math.atan2(node.Y, node.X)
    return x, y

class scanner:
    def __init__(self, model_part, maximum_expected_particle_id,
             outermost_disc_radius, cone_angle, number_of_vanes):
        self.model_part = model_part
        self.n_nodes = maximum_expected_particle_id
        self.disc_radius_2 = outermost_disc_radius ** 2
        self.n_vanes = number_of_vanes
        self.cone_angle = cone_angle
        self.cos_cone_angle_inv = 1. / math.cos(cone_angle)

        zeros   = [0. for i in range(self.n_nodes)]
        self.entered = [0  for i in range(self.n_nodes)]
        self.escaped = [0  for i in range(self.n_nodes)]
        self.ids = []
        self.x0_values = zeros[:]
        self.theta0_values = zeros[:]
        self.t_values  = zeros[:]
        self.x_values  = zeros[:]
        self.y_values  = zeros[:]
        self.z_values  = zeros[:]
        self.vx_values = zeros[:]
        self.vy_values = zeros[:]
        self.vz_values = zeros[:]
        self.wx_values = zeros[:]
        self.wy_values = zeros[:]
        self.wz_values = zeros[:]
        self.D_values  = zeros[:]

    def UpdateData(self, time):

        for node in self.model_part.Nodes:
            i = node.Id

            if node.IsNot(BLOCKED):

                if self.entered[i] == 0:
                    self.entered[i] = 1
                    r, theta = GetPolarCoordinates(node)
                    self.x0_values[i] = r * self.cos_cone_angle_inv
                    self.theta0_values[i] = theta
                    self.D_values[i]  = 2. * node.GetSolutionStepValue(RADIUS, 0)

                r2 = GetPolarRCoordinate2(node)

                if r2 > self.disc_radius_2 and self.escaped[i] == 0:
                    self.ids.append(i)
                    self.escaped[i] = 1
                    self.t_values[i]  = time
                    self.x_values[i]  = node.X
                    self.y_values[i]  = node.Y
                    self.z_values[i]  = node.Z
                    self.vx_values[i] = node.GetSolutionStepValue(VELOCITY_X)
                    self.vy_values[i] = node.GetSolutionStepValue(VELOCITY_Y)
                    self.vz_values[i] = node.GetSolutionStepValue(VELOCITY_Z)
                    self.wx_values[i] = node.GetSolutionStepValue(ANGULAR_VELOCITY_X)
                    self.wy_values[i] = node.GetSolutionStepValue(ANGULAR_VELOCITY_Y)
                    self.wz_values[i] = node.GetSolutionStepValue(ANGULAR_VELOCITY_Z)

    def PrintData(self, in_file_path, out_file_path):
        out_file = open(out_file_path, 'w')

        for node_id in self.ids:
            out_file.write(str( self.t_values[node_id]) + ' ' +
                           str( self.x_values[node_id]) + ' ' +
                           str( self.y_values[node_id]) + ' ' +
                           str( self.z_values[node_id]) + ' ' +
                           str(self.vx_values[node_id]) + ' ' +
                           str(self.vy_values[node_id]) + ' ' +
                           str(self.vz_values[node_id]) + ' ' +
                           str(self.wx_values[node_id]) + ' ' +
                           str(self.wy_values[node_id]) + ' ' +
                           str(self.wz_values[node_id]) + ' ' +
                           str( self.D_values[node_id]) + '\n')
        out_file.close()

        inputs_file = open(in_file_path, 'w')

        for node_id in self.ids:
            inputs_file.write(str(self.x0_values[node_id]) + ' ' +
                              str(self.theta0_values[node_id]) + ' ' +
                              str(self.D_values[node_id]) + '\n')

        inputs_file.close()
