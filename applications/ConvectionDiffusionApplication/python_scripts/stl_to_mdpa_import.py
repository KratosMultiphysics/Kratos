from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import KratosMultiphysics
from KratosMultiphysics.gid_output_process import GiDOutputProcess # Import GiD Output Process

class ModelPartCreator():
    
    def __init__(self, dim_size, input_stl, mdpa_output_file, model_part_name="MainModelPart"):
        
        self._input_stl = input_stl
        self._model_part_name = model_part_name
        self._mdpa_output_file = mdpa_output_file
        self._dimension = dim_size

    def CreateModel(self):
        
        self._initialize_model_part()

        with open (self._input_stl) as read_file:

            node_dict = {}        # dictionary with all vertices. key = (x, y, z); value = node_id
            node_id = 1

            elem_dict = {}
            elem_id = 1

            for row in read_file.readlines():
                row = row.split()

                if (row[0] == "outer"):
                    elem_dict[elem_id] = []
                
                elif (row[0] == "endloop"):
                    elem_id += 1
                
                elif (row[0] == "vertex" and all(self._IsFloat(n) for n in row[1:])):
                    X_coord, Y_coord, Z_coord = [float(coord) for coord in row[1:]]

                    if (X_coord, Y_coord, Z_coord) not in node_dict:
                        node_dict[X_coord, Y_coord, Z_coord] = node_id
                        elem_dict[elem_id].append(node_id)
                        node_id += 1

                    else:
                        for coord, id in node_dict.items():
                            if (coord == (X_coord, Y_coord, Z_coord)):
                                elem_dict[elem_id].append(id)
                                break


        for coord, node_id in node_dict.items():
            
            self.ModelPart.CreateNewNode(node_id, coord[0], coord[1], coord[2])
        
        for id_elem, ids_node in elem_dict.items():
            
            if self._dimension == 2:
                self.ModelPart.CreateNewElement("Element2D3N", id_elem, [ids_node[0], ids_node[1], ids_node[2]], self.ModelPart.GetProperties()[1])
            elif self._dimension == 3:
                Kratos.Logger.PrintWarning("stl_to_mpda_import", "Use it only for 2D models")
            else:
                KratosMultiphysics.Logger.PrintWarning("ModelPartCreator", "Wrong Dimesion Value Prescribed. Choose either (2 or 3).")
        
    def GiDOutput(self, file_name):
        self.gid_output = GiDOutputProcess(self.ModelPart,
                                           file_name,
                                           KratosMultiphysics.Parameters("""
                                  {
                                    "result_file_configuration" : {
                                        "nodal_results"       : []
                                    }
                                  }
                                  """)
                                           )
        self.gid_output.ExecuteInitialize()
        self.gid_output.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteInitializeSolutionStep()
        self.gid_output.PrintOutput()
        self.gid_output.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteFinalize()
    
    def WriteMDPA(self):
        
        KratosMultiphysics.ModelPartIO(self._mdpa_output_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.ModelPart)


    """
    Internal Functions - Not to be accessed by class objects
    """
    
    def _initialize_model_part(self):
        
        current_model = KratosMultiphysics.Model()
        self.ModelPart = current_model.CreateModelPart(self._model_part_name)

        self.ModelPart.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self.ModelPart.SetBufferSize(2)

        self.ModelPart.Nodes.clear()
        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()
    

    def _IsFloat( self, value):
        
        try:
            float( value )
            return True
        except ValueError:
            return False