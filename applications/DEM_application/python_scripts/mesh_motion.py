from __future__ import print_function, absolute_import, division
import math
import os
from KratosMultiphysics import *

def VectSum(v, w):
    return [x + y for x, y in zip(v, w)]

def VectTimes(v, c):
    return [c * x for x in v]

def InnerProd(v, w):
    res = 0.0

    for x, y in zip(v, w):
        res += x * y

    return res

def Norm(v):
    res = 0.0

    for x in v:
        res += x * x

    return math.sqrt(res)

def Normalize(v):

    if all(x == 0.0 for x in v):
        return v

    c = 1 / Norm(v)
    return [x * c for x in v]

def Cross(v, w):
    v_x_w = [v[1] * w[2] - v[2] * w[1],
             v[2] * w[0] - v[0] * w[2],
             v[0] * w[1] - v[1] * w[0]]
    return v_x_w

def RotateRightHandedBasisAroundAxis(e1, e2, axis, ang):
    u = axis[0]
    v = axis[1]
    w = axis[2]
    cang = math.cos(ang)
    sang = math.sin(ang)

    x = e1[0]
    y = e1[1]
    z = e1[2]

    e1[0] = u * (u * x + v * y + w * z) * (1 - cang) + x * cang + (- w * y + v * z) * sang
    e1[1] = v * (u * x + v * y + w * z) * (1 - cang) + y * cang +   (w * x - u * z) * sang
    e1[2] = w * (u * x + v * y + w * z) * (1 - cang) + z * cang + (- v * x + u * y) * sang

    x = e2[0]
    y = e2[1]
    z = e2[2]

    e2[0] = u * (u * x + v * y + w * z) * (1 - cang) + x * cang + (- w * y + v * z) * sang
    e2[1] = v * (u * x + v * y + w * z) * (1 - cang) + y * cang +   (w * x - u * z) * sang
    e2[2] = w * (u * x + v * y + w * z) * (1 - cang) + z * cang + (- v * x + u * y) * sang

    e3 = Cross(e1, e2)

    return [e1, e2, e3]


def MoveAllMeshes(model_part, time):

    if model_part.NumberOfMeshes() > 1:

        for mesh_number in range(1, model_part.NumberOfMeshes()):
            mesh_nodes       = model_part.GetMesh(mesh_number).Nodes
            # print model_part.Properties[mesh_number]
            linear_velocity  = model_part.GetMesh(mesh_number).Properties[0][VELOCITY]
            linear_period    = model_part.GetMesh(mesh_number).Properties[0][VELOCITY_PERIOD]
            angular_velocity = model_part.GetMesh(mesh_number).Properties[0][ANGULAR_VELOCITY]
            angular_period   = model_part.GetMesh(mesh_number).Properties[0][ANGULAR_VELOCITY_PERIOD]
            initial_center   = model_part.GetMesh(mesh_number).Properties[0][ROTATION_CENTER]

            if (linear_period > 0.0):
                linear_omega = 2 * math.pi / linear_period
                center_position = VectSum(initial_center, VectTimes(linear_velocity, math.sin(linear_omega * time) / linear_omega))
            else:
                center_position = initial_center + time * linear_velocity

            if (angular_period > 0.0):
                angular_omega = 2 * math.pi / angular_period
                angle = VectTimes(angular_velocity, math.sin(angular_omega * time) / angular_omega)
            else:
                angle = VectTimes(angular_velocity, time)

            mod_angular_velocity = Norm(angular_velocity)
            relative_position = [0.0, 0.0, 0.0]

            if (mod_angular_velocity > 0.0):
                ang = Norm(angle)
                rotation_axis = Normalize(angular_velocity)
                e1 = [1.0, 0, 0]
                e2 = [0, 1.0, 0]
                [new_axes1, new_axes2, new_axes3] = RotateRightHandedBasisAroundAxis(e1, e2, rotation_axis, ang)

                for node in mesh_nodes:
                    local_X = node.X0 - initial_center[0]
                    local_Y = node.Y0 - initial_center[1]
                    local_Z = node.Z0 - initial_center[2]
                    relative_position[0] = new_axes1[0] * local_X + new_axes2[0] * local_Y + new_axes3[0] * local_Z
                    relative_position[1] = new_axes1[1] * local_X + new_axes2[1] * local_Y + new_axes3[1] * local_Z
                    relative_position[2] = new_axes1[2] * local_X + new_axes2[2] * local_Y + new_axes3[2] * local_Z
                    # NEW POSITION
                    [node.X, node.Y, node.Z] = VectSum(center_position, relative_position)
                    velocity_due_to_rotation = Cross(relative_position, angular_velocity)
                    # NEW VELOCITY                    
                    vel = Vector(3)
                    vel[0] = linear_velocity[0] + velocity_due_to_rotation[0]
                    vel[1] = linear_velocity[1] + velocity_due_to_rotation[1]
                    vel[2] = linear_velocity[2] + velocity_due_to_rotation[2]
                    #The next line only works for Vector(3) or Arrays, not for the lists we are working with here!!
                    node.SetSolutionStepValue(VELOCITY, vel)
                    

