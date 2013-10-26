from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *

CheckForPreviousImport()
import math


def GetNodeBefore(table, prop):
    return table[prop - 1][1]


def GetNodeBegin(table, prop):
    return table[prop - 1][2]


def GetNodeEnd(table, prop):
    return table[prop - 1][3]


def GetNodeAfter(table, prop):
    return table[prop - 1][4]

# select and remove


def DoRemoval(model_part):
    import config
    import full_nodes_table
    print "deactivate_list", config.deactivate_list
    print "list of inlets", config.inlets_1d
    print "list of outlets", config.outlets_1d
    print "full nodes table", full_nodes_table.table
    print model_part

    # mark for deactivation all of the nodes which are not needed
    for prop_id in config.deactivate_list:
        # print 'prop_id', prop_id
        # raw_input()
        for elem in model_part.Elements:
            # print 'properties.id', elem.Properties.Id
            if (elem.Properties.Id == prop_id):
                # print 'hola3'
                elem.SetValue(ERASE_FLAG, True)
                print "..................", elem.Id

                for node in elem.GetNodes():
                    node.SetValue(ERASE_FLAG, True)

    # define a list of nodes which we shall than preserve from removal
    nodes_to_preserve = []
    outlet_nodes = []
    inlet_nodes = []

    # mark for erasal nodes before inlet
    for i in range(0, len(config.inlets_1d)):
        print "nodes before inlet"
        flag_id = config.inlets_1d[i][0]
        prop_id = config.inlets_1d[i][1]
        # print "prop_id", prop_id
        # print GetNodeBefore(full_nodes_table.table, prop_id)
        # raw_input()
        node_before = model_part.Nodes[
            GetNodeBefore(full_nodes_table.table, prop_id)]
        node_begin = model_part.Nodes[
            GetNodeBegin(full_nodes_table.table, prop_id)]
        node_before.SetValue(ERASE_FLAG, True)
        node_begin.SetValue(ERASE_FLAG, True)
        print node_before
        print node_begin
        node_before.SetSolutionStepValue(FLAG_VARIABLE, 0, flag_id)
        print node_begin
        # nodes_to_preserve.append(node_begin)
        nodes_to_preserve.append(node_before)
        inlet_nodes.append(node_before)

    # raw_input()

    # mark for erasal nodes after outlet
    for i in range(0, len(config.outlets_1d)):
        print "nodes after outlet"
        flag_id = config.outlets_1d[i][0]
        prop_id = config.outlets_1d[i][1]
        node_end = model_part.Nodes[
            GetNodeEnd(full_nodes_table.table, prop_id)]
        node_after = model_part.Nodes[
            GetNodeAfter(full_nodes_table.table, prop_id)]
        node_end.SetValue(ERASE_FLAG, True)
        node_after.SetValue(ERASE_FLAG, True)
        print "prop_id = ", prop_id
        print "node end = ", node_end.Id
        print "node after = ", node_after.Id
        node_after.SetSolutionStepValue(FLAG_VARIABLE, 0, flag_id)
        nodes_to_preserve.append(node_after)
        outlet_nodes.append(node_after)
        # raw_input()

    # mark for deactivation the conditions which have all of their nodes
    # marked for erasal
    for cond in model_part.Conditions:
        aaa = 0
        tot_nodes = 0
        for node in cond.GetNodes():
            tot_nodes += 1
            if (node.GetValue(ERASE_FLAG) == True):
                aaa += 1
        if (tot_nodes == aaa and aaa != 0):
            cond.SetValue(ERASE_FLAG, True)
            print "cond", cond
        # print "cond", cond
    print "CONDICIONES A ELIMINAR"
    # raw_input()

    for node in model_part.Nodes:
        flag = node.GetSolutionStepValue(FLAG_VARIABLE)
        if(flag != 0):
            print "id = ", node.Id, " flag = ", flag
    

    # unamrk the node to be preserved
    for node in nodes_to_preserve:
        node.SetValue(ERASE_FLAG, False)
        print "node to preserve", node

    # raw_input()
    print "ANTES DE ELIMINAR NODOS Y CONDICIONES"
    # do delete elements conditions and nodes
    NodeEraseProcess(model_part).Execute()
    ElementEraseProcess(model_part).Execute()
    ConditionEraseProcess(model_part).Execute()

    # raw_input()
    for cond in model_part.Conditions:
        print "cond", cond

    print "DESPUES DE ELIMINAR NODOS Y CONDICIONES"
    # raw_input()
    highest_cond_id = 0
    for cond in model_part.Conditions:
        if (cond.Id > highest_cond_id):
            highest_cond_id = cond.Id
    print highest_cond_id

    # raw_input()
    # add conditions associated to the outlet of the 3d
    id_of_new_property = 0
    for node in inlet_nodes:
        print node
        node.Fix(NODAL_AREA)
        # node.Fix(FLOW)
        # print "flow on outlet node ",node.GetSolutionStepValue(FLOW)
        print "outlet_node =", node.Id
        CreateNewCondition("ArteryOutletCondition", model_part,
                           highest_cond_id + 1, id_of_new_property, node)
        print id_of_new_property
        # print model_part
        print node
        # print highest_cond_id + 1
        highest_cond_id += 1

    id_of_new_property = 0
    for node in outlet_nodes:
        # print node
        print "inlet_node =", node.Id
        node.Fix(FLOW)
        # node.Fix(NODAL_AREA)
        # print model_part
        CreateNewCondition("ArteryInletCondition", model_part,
                           highest_cond_id + 1, id_of_new_property, node)
        # print id_of_new_property
        # print model_part
        # print node
        # print highest_cond_id + 1
        highest_cond_id += 1

    # add conditions associated to the inlet of the 3d
    # for i in range(0, len(config.inlets_1d)):
        # flag_id = config.inlets_1d[i][0]
        # prop_id = config.inlets_1d[i][1]
        # node_before = model_part.Nodes[GetNodeBefore(full_nodes_table.table, prop_id)]
        # node_begin = model_part.Nodes[GetNodeBegin(full_nodes_table.table, prop_id)]
        # CreateNewCondition("Artery1Dto3DCondition", model_part, highest_cond_id + 1, id_of_new_property, node_before,
                           # node_begin)
        # highest_cond_id += 1
        # print node_before
        # print node_begin
        # raw_input()
    for cond in model_part.Conditions:
        print "cond", cond.Id

    print "CREADAS CONDICIONES"
    raw_input()


def ComputePressure(model_part1D):
    initial_pressure = 0
    for node in model_part1D.Nodes:
        # print node.Id
        beta = (
            node.GetSolutionStepValue(YOUNG_MODULUS) * node.GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)) / (
                1.0 - (node.GetSolutionStepValue(POISSON_RATIO) * node.GetSolutionStepValue(POISSON_RATIO)))
        A = node.GetSolutionStepValue(NODAL_AREA)
        A0 = node.GetValue(NODAL_AREA)
        # print A0
        # press=outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
        press = initial_pressure + beta * (math.sqrt(A) - math.sqrt(A0)) / A0
        # press = initial_pressure+beta*(math.sqrt(A/A0)) - beta
        node.SetSolutionStepValue(PRESSURE, 0, press)
