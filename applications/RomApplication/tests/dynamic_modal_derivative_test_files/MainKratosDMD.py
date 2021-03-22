from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.RomApplication as RomApplication
from KratosMultiphysics.RomApplication.structural_mechanics_modal_derivative_analysis import StructuralMechanicsModalDerivativeAnalysis
import json

class StructuralMechanicsModalDerivativeAnalysisMSConstraints(StructuralMechanicsModalDerivativeAnalysis):

    def ModifyInitialGeometry(self):
        # pass
        super().ModifyInitialGeometry()
        model_part = self.model.GetModelPart("Structure")

        constraint_id = 1
        constant = 0.0

        # Top        
        master_nodes_sub_model_part_name = "GENERIC_Points_Beam_Top"
        slave_node_sub_model_part_name = "GENERIC_Points_Top"
        
        master_nodes_sub_model_part = model_part.GetSubModelPart(master_nodes_sub_model_part_name)
        slave_node_sub_model_part = model_part.GetSubModelPart(slave_node_sub_model_part_name)
        
        master_nodes = master_nodes_sub_model_part.Nodes
        slave_nodes = slave_node_sub_model_part.Nodes

        weight = 1.0/float(len(master_nodes))
        for master_node in master_nodes:
            for slave_node in slave_nodes:
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.DISPLACEMENT_X, slave_node, KratosMultiphysics.DISPLACEMENT_X, weight, constant)
                constraint_id += 1
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.DISPLACEMENT_Y, slave_node, KratosMultiphysics.DISPLACEMENT_Y, weight, constant)
                constraint_id += 1
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.DISPLACEMENT_Z, slave_node, KratosMultiphysics.DISPLACEMENT_Z, weight, constant)
                constraint_id += 1

        # Level 1
        master_nodes_sub_model_part_name = "GENERIC_Points_Surface_Mid"
        slave_node_sub_model_part_name = "GENERIC_Points_Truss_Mid"
        
        master_nodes_sub_model_part = model_part.GetSubModelPart(master_nodes_sub_model_part_name)
        slave_node_sub_model_part = model_part.GetSubModelPart(slave_node_sub_model_part_name)
        
        master_nodes = master_nodes_sub_model_part.Nodes
        slave_nodes = slave_node_sub_model_part.Nodes

        weight = 1.0/float(len(master_nodes))
        for master_node in master_nodes:
            for slave_node in slave_nodes:
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.DISPLACEMENT_X, slave_node, KratosMultiphysics.DISPLACEMENT_X, weight, constant)
                constraint_id += 1
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.DISPLACEMENT_Y, slave_node, KratosMultiphysics.DISPLACEMENT_Y, weight, constant)
                constraint_id += 1
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.DISPLACEMENT_Z, slave_node, KratosMultiphysics.DISPLACEMENT_Z, weight, constant)
                constraint_id += 1

        # Level 2
        master_nodes_sub_model_part_name = "GENERIC_Points_Shell_Mid"
        slave_node_sub_model_part_name = "GENERIC_Points_Beam_Mid"
        
        master_nodes_sub_model_part = model_part.GetSubModelPart(master_nodes_sub_model_part_name)
        slave_node_sub_model_part = model_part.GetSubModelPart(slave_node_sub_model_part_name)
        
        master_nodes = master_nodes_sub_model_part.Nodes
        slave_nodes = slave_node_sub_model_part.Nodes

        weight = 1.0/float(len(master_nodes))
        for master_node in master_nodes:
            for slave_node in slave_nodes:
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.DISPLACEMENT_X, slave_node, KratosMultiphysics.DISPLACEMENT_X, weight, constant)
                constraint_id += 1
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.DISPLACEMENT_Y, slave_node, KratosMultiphysics.DISPLACEMENT_Y, weight, constant)
                constraint_id += 1
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.DISPLACEMENT_Z, slave_node, KratosMultiphysics.DISPLACEMENT_Z, weight, constant)
                constraint_id += 1
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.ROTATION_X, slave_node, KratosMultiphysics.ROTATION_X, weight, constant)
                constraint_id += 1
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.ROTATION_Y, slave_node, KratosMultiphysics.ROTATION_Y, weight, constant)
                constraint_id += 1
                model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, master_node, KratosMultiphysics.ROTATION_Z, slave_node, KratosMultiphysics.ROTATION_Z, weight, constant)
                constraint_id += 1

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsModalDerivativeAnalysisMSConstraints(model,parameters)
    simulation.Run()
