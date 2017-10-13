from __future__ import print_function, absolute_import, division
import math
import os
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

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

def ApplyEmbeddedBCsToFluid(model_part):

    for node in model_part.Nodes:
        old_dist = node.GetSolutionStepValue(DISTANCE, 1)
        dist     = node.GetSolutionStepValue(DISTANCE)

        if (dist < 0.0):
            node.Fix(PRESSURE)
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.Fix(VELOCITY_Z)
            node.SetSolutionStepValue(PRESSURE, 0.0)
            node.SetSolutionStepValue(VELOCITY_X, 0.0)
            node.SetSolutionStepValue(VELOCITY_Y, 0.0)
            node.SetSolutionStepValue(VELOCITY_Z, 0.0)

        elif (old_dist < 0.0):
            node.Free(PRESSURE)
            node.Free(VELOCITY_X)
            node.Free(VELOCITY_Y)
            node.Free(VELOCITY_Z)
            
    if model_part.NumberOfMeshes() > 1:

        for mesh_number in range(2, model_part.NumberOfMeshes()):
            mesh_nodes = model_part.GetMesh(mesh_number).Nodes
            # print model_part.Properties[mesh_number]

            # INLETS
            if (model_part.Properties[mesh_number][IMPOSED_VELOCITY_X] == 1 or model_part.Properties[mesh_number][IMPOSED_VELOCITY_Y] == 1 or model_part.Properties[mesh_number][IMPOSED_VELOCITY_Z] == 1):

                for node in mesh_nodes:
                    dist = node.GetSolutionStepValue(DISTANCE)

                    if (dist > 0.0):
                        node.Free(PRESSURE)
                        if (model_part.Properties[mesh_number][IMPOSED_VELOCITY_X]):
                            node.Fix(VELOCITY_X)
                            node.SetSolutionStepValue(VELOCITY_X, model_part.Properties[mesh_number][IMPOSED_VELOCITY_X_VALUE])
                        if (model_part.Properties[mesh_number][IMPOSED_VELOCITY_Y]):
                            node.Fix(VELOCITY_Y)
                            node.SetSolutionStepValue(VELOCITY_Y, model_part.Properties[mesh_number][IMPOSED_VELOCITY_Y_VALUE])
                        if (model_part.Properties[mesh_number][IMPOSED_VELOCITY_Z]):
                            node.Fix(VELOCITY_Z)
                            node.SetSolutionStepValue(VELOCITY_Z, model_part.Properties[mesh_number][IMPOSED_VELOCITY_Z_VALUE])

            # OUTLETS
            if (model_part.Properties[mesh_number][IMPOSED_PRESSURE] == 1):
                # here I assume all nodes of this outlet have the same body force and density!!

                for node in mesh_nodes:
                    bf = node.GetSolutionStepValue(BODY_FORCE)
                    mod_bf = Norm(bf)
                    normalized_bf = Normalize(bf)
                    outlet_density = node.GetSolutionStepValue(DENSITY)
                    break

                for node in mesh_nodes:
                    #coord = Vector(node)
                    height = - InnerProd([node.X, node.Y, node.Z], normalized_bf)
                    maxheight = -1.0e90

                    if (height > maxheight):
                        maxheight = height
                        highest_node = node

                base_pressure = model_part.Properties[mesh_number][PRESSURE]

                for node in mesh_nodes:
                    dist = node.GetSolutionStepValue(DISTANCE)
                    distance_to_highest = maxheight - InnerProd(VectTimes([node.X, node.Y, node.Z], -1), normalized_bf)

                    if (dist > 0.0):
                        node.Fix(PRESSURE)
                        actual_pressure = base_pressure + outlet_density * mod_bf * distance_to_highest
                        node.SetSolutionStepValue(PRESSURE, actual_pressure)
                        node.Free(VELOCITY_X)
                        node.Free(VELOCITY_Y)
                        node.Free(VELOCITY_Z)

def ApplyEmbeddedBCsToBalls(model_part, DEMParameters):

    for node in model_part.Nodes:
        dist = node.GetSolutionStepValue(DISTANCE)

        if (dist < 0.0):

            if node.Is(BLOCKED):
                node.Set(ACTIVE, True)

            elif (DEMParameters.RemoveBallsInEmbeddedOption):
                #print("--------------- One particle was erased because it entered an embedded structure -------------------")
                node.Set(TO_ERASE,True)


class SearchEmbeddedDEMTools:

    def __init__(self):
        self.search_tools = DemSearchUtilities(OMP_DEMSearch())

    def SearchNodeNeighboursDistances(self, model_part, dem_model_part, search_radius):
        self.search_tools.SearchNodeNeighboursDistances(model_part, dem_model_part, search_radius, DISTANCE)

    def CalculateElementNeighbourDistances(self, model_part, intersecting_surface_semi_thickness):

        for node in model_part.Nodes:
            distance = node.GetSolutionStepValue(DISTANCE) - intersecting_surface_semi_thickness
            node.SetSolutionStepValue(DISTANCE, 0, distance)

            if (distance < 0):
                node.Fix(PRESSURE)
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Y)
                node.Fix(VELOCITY_Z)

        for element in model_part.Elements:
            negative = 0
            positive = 0

            for node in element.GetNodes():
                d = node.GetSolutionStepValue(DISTANCE)

                if (d >= 0.0):
                    positive = positive + 1

                else:
                    negative = negative + 1

            if ((negative > 0) and (positive > 0)):
                tmp = Vector(4)
                i = 0
                element.SetValue(SPLIT_ELEMENT, True)

                for node in element.GetNodes():
                    d = node.GetSolutionStepValue(DISTANCE)
                    tmp[i] = d
                    i = i + 1

                element.SetValue(ELEMENTAL_DISTANCES, tmp)
