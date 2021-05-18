import numpy
from stl import mesh #this requires numpy-stl
import KratosMultiphysics

# Using an existing stl file:
input_name = "Engine.stl"
mesh = mesh.Mesh.from_multi_file(input_name)
#print(mesh.Mesh)
my_mp = KratosMultiphysics.ModelPart("skin")
prop = my_mp.Properties[0]

node_id = 1
elem_id = 1
for m in mesh:
    print(m)

    for vertex in m.points:
        n1 = my_mp.CreateNewNode(node_id, float(vertex[0]), float(vertex[1]), float(vertex[2]))
        node_id+=1
        n2 = my_mp.CreateNewNode(node_id, float(vertex[3]), float(vertex[4]), float(vertex[5]))
        node_id+=1
        n3 = my_mp.CreateNewNode(node_id, float(vertex[6]), float(vertex[7]), float(vertex[8]))
        node_id+=1

        el = my_mp.CreateNewElement("Element3D3N",elem_id,  [n1.Id, n2.Id, n3.Id], prop)
        elem_id += 1

print(my_mp)

from gid_output_process import GiDOutputProcess
gid_output = GiDOutputProcess(my_mp,
                              input_name,
                              KratosMultiphysics.Parameters("""
                                  {
                                    "result_file_configuration" : {
                                        "nodal_results"       : []
                                    }
                                  }
                                  """)
                              )
gid_output.ExecuteInitialize()
gid_output.ExecuteBeforeSolutionLoop()

gid_output.ExecuteInitializeSolutionStep()
gid_output.PrintOutput()
gid_output.ExecuteFinalizeSolutionStep()
gid_output.ExecuteFinalize()
