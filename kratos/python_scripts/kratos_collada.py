# importing the Kratos Library
from KratosMultiphysics import *

#importing the collada library
import collada
import numpy

class ColladaImporter:
    def __init__(self,model_part, collada_filename, import_normals=False):
        self.model_part = model_part
        self.collada_filename = str(collada_filename)
        self.import_normals = import_normals

        #open and import the collada collada_file
        print("attempting to open ", self.collada_filename)
        f = open(self.collada_filename,'rb')

        self.collada_mesh = collada.Collada(f)
        f.close()

    def GenerateModelPartElements(self, element_name):
        #kratos auxiliaries
        nodeid = 1
        elid = 1
        prop = self.model_part.Properties[0]

        ###here loop over the skins in needed
        print("controller = ",self.collada_mesh.controllers)
        for skin in self.collada_mesh.controllers:
            print("skin joint_matrices",skin.joint_matrices)
            print(" ")

        ##following block reads in a collada static file and generates the kratos database
        scene = self.collada_mesh.scene
        for mesh in scene.objects('geometry'):
            print("mesh = ",mesh)
            #first of all generate all of the kratos nodes for this mesh
            primitives_list = list(mesh.primitives())
            #print(primitives_list)
            for primitive in primitives_list:
                print("primitive type = ", type(primitive) )
                if( type(primitive) == collada.triangleset.BoundTriangleSet):
                    triset = primitive
                elif( type(primitive) == collada.polylist.BoundPolylist):
                    triset = primitive.triangleset()
                else:
                    triset = None


                if triset != None:

                    #here we loop over all of the triangles
                    for tri_indices,normal_indices in zip(triset.vertex_index,triset.normal_index):
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
                        #print("just before creating the element")
                        self.model_part.CreateNewElement(element_name, elid,  elemental_ids, prop)
                        elid += 1

                        if self.import_normals == True:
                            ##here we assign the normals
                            n0a = triset.normal[normal_indices][0]
                            n0 = Vector(3)
                            n0[0] = float(n0a[0])
                            n0[1] = float(n0a[1])
                            n0[2] = float(n0a[2])

                            n1a = triset.normal[normal_indices][1]
                            n1 = Vector(3)
                            n1[0] = float(n1a[0]); n1[1] = float(n1a[1]); n1[2] = float(n1a[2]);
                            n2a = triset.normal[normal_indices][2]
                            n2 = Vector(3)
                            n2[0] = float(n2a[0]); n2[1] = float(n2a[1]); n2[2] = float(n2a[2]);

                            kratos_node0.SetSolutionStepValue(NORMAL,0,n0)
                            kratos_node1.SetSolutionStepValue(NORMAL,0,n1)
                            kratos_node2.SetSolutionStepValue(NORMAL,0,n2)













