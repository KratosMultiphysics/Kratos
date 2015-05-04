from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# importing the Kratos Library
#from KratosMultiphysics import *
#from KratosMultiphysics.IncompressibleFluidApplication import *
#from KratosMultiphysics.FluidDynamicsApplication import *
#CheckForPreviousImport()

#importing the collada library
import collada
import numpy

class ColladaImporter:
    def __init__(self,model_part, collada_filename):
        self.model_part = model_part
        self.collada_filename = str(collada_filename)
        
        #open and import the collada collada_file
        print("attempting to open ", self.collada_filename)
        f = open(self.collada_filename,'rb')
        print("line 21")
        self.collada_mesh = collada.Collada(f)
        f.close()
    
    def GenerateModelPartConditions(self, condition_name):
        #loop over all the triangles and create new nodes assigning an Id in the order in which they appear
        #at the same time also generate the Conditions
        nodeid = 1
        elid = 1
        prop = self.model_part.Properties[0]
        
        for itmesh in range(0, len(self.collada_mesh.geometries)):
            mesh = self.collada_mesh.geometries[itmesh]
            print("mesh = ",mesh)
            #first of all generate all of the kratos nodes for this mesh
            primitives_list = list(mesh.primitives)
            for triset in primitives_list:
                if( type(triset) == collada.triangleset.TriangleSet): 
                    print("triset = ", type(triset) )
                
                    #here we loop over all of the triangles
                    for tri_indices in triset.vertex_index:
                        node0_coords = triset.vertex[tri_indices][0]
                        node1_coords = triset.vertex[tri_indices][1]
                        node2_coords = triset.vertex[tri_indices][2]
                        #print("node0", node0_coords)
                        #print("node1", node1_coords)
                        #print("node2", node2_coords)
                        
                        ##here we shall generate the kratos node for each of those
                        elemental_ids = []
                        kratos_node0 = self.model_part.CreateNewNode( nodeid, float(node0_coords[0]), float(node0_coords[1]), float(node0_coords[2]))
                        elemental_ids.append(nodeid)
                        nodeid += 1
                        
                        kratos_node1 = self.model_part.CreateNewNode( nodeid, float(node1_coords[0]), float(node1_coords[1]), float(node1_coords[2]))
                        elemental_ids.append(nodeid)
                        nodeid += 1
                        
                        kratos_node2 = self.model_part.CreateNewNode( nodeid, float(node2_coords[0]), float(node2_coords[1]), float(node2_coords[2]))
                        elemental_ids.append(nodeid)
                        nodeid += 1
                        
                        ##here generate the conditions
                        self.model_part.CreateNewCondition(condition_name, elid,  elemental_ids, prop)
                        elid += 1
                        
                        ##we shall also keep this in a tuple for future use, maintaining the order of creation
                    
                    
            
            
            #then generate all of the triangles
            
            
        
        
        
        
    
        