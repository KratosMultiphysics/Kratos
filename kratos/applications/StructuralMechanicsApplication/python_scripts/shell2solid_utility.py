from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from numpy import *

import sys
import os
import os.path
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import * 
CheckForPreviousImport()

def CalculateNormalElement(elem):
    nodes =[]
    for node in elem.GetNodes():
        nodes.append(node)
        
    x = [nodes[1].X - nodes[0].X , nodes[1].Y - nodes[0].Y, nodes[1].Z - nodes[0].Z]
    y = [nodes[2].X - nodes[0].X , nodes[2].Y - nodes[0].Y, nodes[2].Z - nodes[0].Z]
    
    return cross(x, y)


def Shell2SolidUtility(working_folder, output_name, thickness, layer_number,prop_vector):
    
    #### NOTE: It always extrudes in the normal direction!!!!!!
    os.chdir(working_folder)
    # defining the model parts
    model_part = ModelPart("Temp_Shell")
    
    # add displacements
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    # add specific variables for the problem conditions
    model_part.AddNodalSolutionStepVariable(POINT_LOAD)
    
    # reading the shell mesh
    input_file_name = "shell_mesh"
    model_part_io=ModelPartIO(input_file_name)
    model_part_io.ReadModelPart(model_part)
    
    if os.path.isfile(output_name+".mdpa"):
        os.remove(output_name+".mdpa")
        
    file = open(output_name+".mdpa", "w")
    # Write PROPERTIES
    file.write("Begin ModelPartData\n")
    file.write("//VARIABLE_NAME value\n")
    file.write("End ModelPartData\n")
    
    file.write("\n")
    
    file.write("Begin Properties  1\n")
    file.write("CONSTITUTIVE_LAW_NAME  Elastic3DLaw\n")
    file.write("YOUNG_MODULUS"+"          "+prop_vector[0]+"\n")
    file.write("POISSON_RATIO"+"             "+prop_vector[1]+"\n")
    file.write("DENSITY"+"                         "+prop_vector[2]+"\n")
    file.write("NINT_TRANS"+"                   "+prop_vector[3]+"\n")
    file.write("EAS_IMP"+"                         "+prop_vector[4]+"\n")
    file.write("QUAD_ON"+"                      "+prop_vector[5]+"\n")
    file.write("End Properties\n")
    
    file.write("\n")
    
    number_nodes = 0
    for node in model_part.Nodes:
        number_nodes += 1
    
    number_elements = 0
    for element in model_part.Elements:
        number_elements += 1
    
    ## Write NODES
    file.write("Begin Nodes\n")
    
    for node in model_part.Nodes:
        # Look for the normal
        for element in model_part.Elements:
            for elem_nodes in element.GetNodes():
                if node.Id == elem_nodes.Id:
                    normal_elem = CalculateNormalElement(element)
                    
        for layer in range(layer_number + 1):
            aux = layer/layer_number * thickness
            file.write(str(node.Id + layer * number_nodes)+"\t"+str(node.X + aux * normal_elem[0])+"\t"+str(node.Y + aux * normal_elem[1])+"\t"+str(node.Z + aux * normal_elem[2])+"\n")
        
    file.write("End Nodes\n")

    file.write("\n")

    ## Write ELEMENT
    file.write("Begin Elements SprismElement3D6N\n")
    
    for element in model_part.Elements:
        nodes = element.GetNodes()
        for layer in range(layer_number):
            file.write(str(element.Id + layer * number_elements)+"\t 1 \t"+str(nodes[0].Id + (layer + 1) * number_nodes)+"\t"+str(nodes[1].Id + (layer + 1) * number_nodes)+"\t"+str(nodes[2].Id + (layer + 1) * number_nodes)+"\t"+str(nodes[0].Id + layer * number_nodes)+"\t"+str(nodes[1].Id+ layer * number_nodes)+"\t"+str(nodes[2].Id + layer * number_nodes)+"\n")
            
    file.write("End Elements\n")
    
    file.write("\n")
    
    ## Write BC
    # Write DISP_X
    file.write("Begin NodalData DISPLACEMENT_X\n")
    for node in model_part.Nodes:
        if(node.IsFixed(DISPLACEMENT_X)):
            for layer in range(layer_number):
                file.write(str(node.Id + layer * number_nodes)+"\t1\t0.00000\n")
    file.write("End NodalData\n")
    
    file.write("\n")
    
    # Write DISP_Y
    file.write("Begin NodalData DISPLACEMENT_Y\n")
    for node in model_part.Nodes:
        if(node.IsFixed(DISPLACEMENT_Y)):
            for layer in range(layer_number):
                file.write(str(node.Id + layer * number_nodes)+"\t1\t0.00000\n")
    file.write("End NodalData\n")
    
    file.write("\n")
    
    # Write DISP_Z
    file.write("Begin NodalData DISPLACEMENT_Z\n")
    for node in model_part.Nodes:
        if(node.IsFixed(DISPLACEMENT_Z)):
            for layer in range(layer_number):
                file.write(str(node.Id + layer * number_nodes)+"\t1\t 0.00000\n")
    file.write("End NodalData\n")
    
    file.write("\n")
    
    ## Write POINT_LOAD
    file.write("Begin Conditions PointLoadCondition3D1N\n")
    counter = 0
    for node in model_part.Nodes:
        load = node.GetSolutionStepValue(POINT_LOAD, 0)
        if ((abs(load[0]) > 0.0) | (abs(load[1]) > 0.0) | (abs(load[2]) > 0.0)):
            for layer in range(layer_number):
                counter += 1
                file.write(str(counter)+"\t1\t"+str(node.Id + layer * number_nodes)+"\n")
    file.write("End Conditions\n")
    
    file.write("\n")
    
    # Write LOAD_X
    file.write("Begin NodalData POINT_LOAD_X\n")
    for node in model_part.Nodes:
        load = node.GetSolutionStepValue(POINT_LOAD, 0)
        if (abs(load[0]) > 0):
            for layer in range(layer_number):
                 file.write(str(node.Id + layer * number_nodes)+"\t1\t"+str(load[0])+"\n")
    file.write("End NodalData\n")
    
    file.write("\n")
    
    # Write LOAD_Y
    file.write("Begin NodalData POINT_LOAD_Y\n")
    for node in model_part.Nodes:
        load = node.GetSolutionStepValue(POINT_LOAD, 0)
        if (abs(load[1]) > 0):
            for layer in range(layer_number):
                 file.write(str(node.Id + layer * number_nodes)+"\t1\t"+str(load[1])+"\n")
    file.write("End NodalData\n")
    
    file.write("\n")
    
    # Write LOAD_Z
    file.write("Begin NodalData POINT_LOAD_Z\n")
    for node in model_part.Nodes:
        load = node.GetSolutionStepValue(POINT_LOAD, 0)
        if (abs(load[2]) > 0):
            for layer in range(layer_number):
                 file.write(str(node.Id + layer * number_nodes)+"\t1\t"+str(load[2])+"\n")
    file.write("End NodalData\n")
    
    file.write("\n")
    
    file.close()
    


