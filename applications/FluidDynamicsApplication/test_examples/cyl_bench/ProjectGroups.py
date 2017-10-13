from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from Kratos import *


# This part to be configured for each problem


project_groups_elements_ids = {}

# This part is always the same for all problems

model_groups = {}


def InitializeGroups(model_part):
    for group_name in project_groups_names:
        mesh = Mesh()

        for node_id in project_groups_nodes_ids[group_name]:
            mesh.Nodes[node_id] = model_part.Nodes[node_id]

        for element_id in project_groups_elements_ids[group_name]:
            mesh.Elements[element_id] = model_part.Elements[element_id]

        for condition_id in project_groups_conditions_ids[group_name]:
            mesh.Conditions[condition_id] = model_part.Conditions[condition_id]

        print(group_name, "mesh:", mesh)

        model_groups[group_name] = mesh
