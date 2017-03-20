from __future__ import print_function, absolute_import, division
import math
import os
from KratosMultiphysics import *

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
    
    cang         = math.cos(ang)
    sang         = math.sin(ang)
    
    new_axes1    = Vector(3)
    new_axes1[0] = axis[0] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[0] * cang + (- axis[2] * e1[1] + axis[1] * e1[2]) * sang
    new_axes1[1] = axis[1] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[1] * cang +   (axis[2] * e1[0] - axis[0] * e1[2]) * sang
    new_axes1[2] = axis[2] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[2] * cang + (- axis[1] * e1[0] + axis[0] * e1[1]) * sang
    
    new_axes2    = Vector(3)   
    new_axes2[0] = axis[0] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[0] * cang + (- axis[2] * e2[1] + axis[1] * e2[2]) * sang
    new_axes2[1] = axis[1] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[1] * cang +   (axis[2] * e2[0] - axis[0] * e2[2]) * sang
    new_axes2[2] = axis[2] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[2] * cang + (- axis[1] * e2[0] + axis[0] * e2[1]) * sang    

    new_axes3    = Vector(3)
    new_axes3    = Cross(new_axes1, new_axes2)

    return [new_axes1, new_axes2, new_axes3]

def MoveAllMeshes(model_part, time):

    if model_part.NumberOfMeshes() > 0:

        for mesh_number in range(1, model_part.NumberOfMeshes()):
            mesh_nodes         = model_part.GetMesh(mesh_number).Nodes
            linear_velocity    = model_part.GetMesh(mesh_number)[VELOCITY]
            linear_period      = model_part.GetMesh(mesh_number)[VELOCITY_PERIOD]
            angular_velocity   = model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY]
            angular_period     = model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY_PERIOD]
            initial_center     = model_part.GetMesh(mesh_number)[ROTATION_CENTER]
            
            center_position    = Vector(3)
            center_position[0] = 0.0
            center_position[1] = 0.0
            center_position[2] = 0.0 
            
            if (linear_period > 0.0):
                linear_omega = 2 * math.pi / linear_period
                inv_linear_omega = 1/linear_omega
                center_position = initial_center + linear_velocity * math.sin(linear_omega * time)* inv_linear_omega
                linear_velocity_changed = linear_velocity * math.cos(linear_omega * time)
                
            else:
                center_position = initial_center + time * linear_velocity
                linear_velocity_changed = linear_velocity

            if (angular_period > 0.0):
                angular_omega = 2 * math.pi / angular_period
                inv_angular_omega = 1/angular_omega
                angle = angular_velocity * math.sin(angular_omega * time) * inv_angular_omega
                angular_velocity_changed = angular_velocity * math.cos(angular_omega * time)
                
            else:
                angle = angular_velocity * time
                angular_velocity_changed = angular_velocity
                            
            mod_angular_velocity = Norm(angular_velocity)
            relative_position    = Vector(3)
            relative_position[0] = 0.0
            relative_position[1] = 0.0
            relative_position[2] = 0.0                                    
            
            if (mod_angular_velocity > 0.0):
                ang = Norm(angle)
                rotation_axis = Normalize(angular_velocity)
                e1 = Vector(3)
                e1[0] = 1.0
                e1[1] = 0.0
                e1[2] = 0.0
                
                e2 = Vector(3)
                e2[0] = 0.0
                e2[1] = 1.0
                e2[2] = 0.0
                
                [new_axes1, new_axes2, new_axes3] = RotateRightHandedBasisAroundAxis(e1, e2, rotation_axis, ang)  
                
            else:
                new_axes1    = Vector(3)
                new_axes1[0] = 1.0
                new_axes1[1] = 0.0
                new_axes1[2] = 0.0
                
                new_axes2    = Vector(3)
                new_axes2[0] = 0.0
                new_axes2[1] = 1.0
                new_axes2[2] = 0.0
                
                new_axes3    = Vector(3)
                new_axes3[0] = 0.0
                new_axes3[1] = 0.0
                new_axes3[2] = 1.0
                        
            if (mod_angular_velocity>0.0 or Norm(linear_velocity)>0.0):

                for node in mesh_nodes:
                    local_X = node.X0 - initial_center[0]
                    local_Y = node.Y0 - initial_center[1]
                    local_Z = node.Z0 - initial_center[2]
                    
                    relative_position[0] = new_axes1[0] * local_X + new_axes2[0] * local_Y + new_axes3[0] * local_Z
                    relative_position[1] = new_axes1[1] * local_X + new_axes2[1] * local_Y + new_axes3[1] * local_Z
                    relative_position[2] = new_axes1[2] * local_X + new_axes2[2] * local_Y + new_axes3[2] * local_Z
                    
                    # NEW POSITION
                    [node.X, node.Y, node.Z] = center_position + relative_position
                    
                    displacement = Vector(3)
                    displacement[0] = node.X - node.X0
                    displacement[1] = node.Y - node.Y0
                    displacement[2] = node.Z - node.Z0
                    
                    velocity_due_to_rotation = Cross(angular_velocity_changed , relative_position)
                    
                    # NEW VELOCITY                    
                    vel = Vector(3)
                    vel[0] = linear_velocity_changed[0] + velocity_due_to_rotation[0]
                    vel[1] = linear_velocity_changed[1] + velocity_due_to_rotation[1]
                    vel[2] = linear_velocity_changed[2] + velocity_due_to_rotation[2]
                    
                    #The next line only works for Vector(3) or Arrays, not for the lists we are working with here!!
                    
                    #update VELOCITY
                    node.SetSolutionStepValue(VELOCITY, vel)
                    #update DISPLACEMENT
                    node.SetSolutionStepValue(DISPLACEMENT, displacement)
