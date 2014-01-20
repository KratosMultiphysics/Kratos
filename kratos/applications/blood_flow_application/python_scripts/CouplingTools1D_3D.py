from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
CheckForPreviousImport()

import math
import time


def InitializeInletNodes(model_part_1d, model_part_3d):
    inlet_nodes_1d = []
    inlet_nodes_3d = []
    inlet_area_3d = []

    for node in model_part_1d.Nodes:
        if(node.GetSolutionStepValue(FLAG_VARIABLE) == 1):
            inlet_nodes_1d.append(node)

    for node in model_part_3d.Nodes:
        if(node.GetSolutionStepValue(FLAG_VARIABLE) == 1):
            inlet_nodes_3d.append(node)

    print(inlet_nodes_3d)
    # measure 3D area

    BodyNormalCalculationUtils().CalculateBodyNormals(model_part_3d, 3)
    area3d = 0.0
    for node in inlet_nodes_3d:
        n = node.GetSolutionStepValue(NORMAL)
        a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
        area3d += a

    # area3d = 0.0
    # for cond in model_part_3d.Conditions:
        # n = node.GetSolutionStepValue(NORMAL)
        # a = math.sqrt(n[0]**2 +  n[1]**2 + n[2]**2 )
        # print a
        # for node in cond.GetNodes():
            # if(node.GetSolutionStepValue(FLAG_VARIABLE) == True):
            # area3d += a*0.3333333333333333333333333333
    print("************************,area 3D", area3d)

    # fix conditions as needed
    for node in inlet_nodes_1d:
        node.SetValue(NODAL_AREA, area3d)
        node.SetSolutionStepValue(NODAL_AREA, 0, area3d)
        # node.Fix(PRESSURE)
        node.Fix(NODAL_AREA)
        node.Fix(FLOW)
        # node.SetValue(NODAL_AREA,area3d)
        # node.SetSolutionStepValue(NODAL_AREA,area3d)

    for node in inlet_nodes_3d:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
        node.Fix(VELOCITY_Z)

    return [inlet_nodes_1d, inlet_nodes_3d, area3d]


def InitializeOutletNodes(model_part_1d, model_part_3d):
    outlet_nodes_1d = []
    outlet_nodes_3d = []

    for node in model_part_1d.Nodes:
        if(node.GetSolutionStepValue(FLAG_VARIABLE) == 2):
            outlet_nodes_1d.append(node)
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.Fix(VELOCITY_Z)

    for node in model_part_3d.Nodes:
        if(node.GetSolutionStepValue(FLAG_VARIABLE) == 2):
            outlet_nodes_3d.append(node)
            node.Fix(PRESSURE)

    return [outlet_nodes_1d, outlet_nodes_3d]


def Transfer1D_to_3D(model_part_1d, model_part_3d, inlet_nodes_1d, inlet_nodes_3d, outlet_nodes_1d, outlet_nodes_3d, area3d):

    vel1d = inlet_nodes_1d[0].GetSolutionStepValue(FLOW) / area3d
    print("vel1d = ", vel1d)

    for node in inlet_nodes_3d:
            # aux = node.GetSolutionStepValue(NORMAL);
            # A = math.sqrt(aux[0]*aux[0] + aux[1]*aux[1] + aux[2]*aux[2])
            # if(A == 0.0):
            # print "node id = ",node.Id
            # aux *= vel1d/A
        node.SetSolutionStepValue(VELOCITY_X, 0, vel1d)

    # assign pressure to outlet
    # for i in range(0,len(outlet_nodes_1d) ):
        # press = outlet_nodes_1d[i].GetSolutionStepValue(PRESSURE)
    beta = outlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * outlet_nodes_1d[
        0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)
    A = outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
    A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
    # print "in Transfer1D_to_3D",A0
    press = math.sqrt(A / A0) * beta - beta
    # Edu
    # press = beta/A0*(math.sqrt(A) - math.sqrt(A0) )
    print("in Transfer1D_to_3D", press)

    for node in outlet_nodes_3d:
        node.SetSolutionStepValue(PRESSURE, 0, press)


def Transfer3D_to_1D(model_part_1d, model_part_3d, inlet_nodes_1d, inlet_nodes_3d, outlet_nodes_1d, outlet_nodes_3d, area3d):
    press_3d = 0.0
    counter = 0.0
    for node in inlet_nodes_3d:
        press_3d += node.GetSolutionStepValue(PRESSURE)
        counter += 1.0
    avg_press = press_3d / counter

    # TODO: make it to read from the input
    print("inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS)", inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS))
    print("inlet_nodes_1d[0].GetSolutionStepValue(THICKNESS)", inlet_nodes_1d[0].GetSolutionStepValue(THICKNESS))
    print(inlet_nodes_1d[0].Id)
    beta = inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * inlet_nodes_1d[
        0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)
    A = area3d * (avg_press / beta + 1) ** 2

    # Edu
    # A = (avg_press*A0/beta + math.sqrt(A0) )**2
    print("in Transfer3D_to_1D", A)

    inlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA, 0, A)

    # print "Area on 3D inlet ",A

    # TODO: transform pressure to area

    # assign flow to the outlet
    flow = 0.0
    for node in outlet_nodes_3d:
        normal = node.GetSolutionStepValue(NORMAL)
        vel = node.GetSolutionStepValue(VELOCITY)
        flow += normal[0] * vel[0] + normal[1] * vel[1] + normal[2] * vel[2]

    # print "flow = ",flow
    # Edu
    print("sin Riemman", flow)
    # flow=2*flow-outlet_nodes_1d[0].GetSolutionStepValue(FLOW, 1)
    # print "con Riemman",flow

    outlet_nodes_1d[0].SetSolutionStepValue(FLOW, 0, flow)
