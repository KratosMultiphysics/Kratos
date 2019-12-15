from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy
from stl import mesh #this requires numpy-stl
import write_utility as wf
from KratosMultiphysics import *
from gid_output_process import GiDOutputProcess
        

class MDPAWriter():

    def __init__(self, inputSTL):
      self.input_name = inputSTL
      self.my_mesh = mesh.Mesh.from_file(self.input_name)

    def CreateModel(self, modelName):
      CurrentModel = Model()
      self.myModel = CurrentModel.CreateModelPart(modelName)
      self.prop = self.myModel.Properties[0]

      node_id = 1
      elem_id = 1
      for vertex in self.my_mesh:
          print(vertex)
          n1 = self.myModel.CreateNewNode(node_id, float(vertex[0]), float(vertex[1]), float(vertex[2]))
          node_id+=1
          n2 = self.myModel.CreateNewNode(node_id, float(vertex[3]), float(vertex[4]), float(vertex[5]))
          node_id+=1
          n3 = self.myModel.CreateNewNode(node_id, float(vertex[6]), float(vertex[7]), float(vertex[8]))
          node_id+=1

          el = self.myModel.CreateNewElement("Element2D3N",elem_id,  [n1.Id, n2.Id, n3.Id], self.prop)
          elem_id += 1
      
      print(self.myModel)
    
    def GiDOutput(self):
      self.gid_output = GiDOutputProcess(self.myModel,
                              self.input_name,
                              Parameters("""
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

    def Write_MDPA(self, filename, element_name):
      mdpa_file = open(filename + ".mdpa", 'w')
      wf.PrintModelData(mdpa_file)
      wf.PrintProperties(mdpa_file)
      wf.PrintNodes(self.myModel.Nodes, mdpa_file)
      wf.PrintElements(element_name, self.myModel.Elements, mdpa_file)
      #wf.PrintConditions(condition_name, self.myModel.Conditions, mdpa_file)
      mdpa_file.close()


if __name__ == "__main__":
    inputSTL = "mycircle2.stl" # Enter input STL file name
    modelName = "solid" # Enter Model Name
    filename = "solid_circle"
    element_name = "Element2D3N"

    Writer = MDPAWriter(inputSTL)
    Writer.CreateModel(modelName)
    Writer.Write_MDPA(filename, element_name)
    #Writer.GiDOutput()

#if __name__ == "__create_mdpa__":
#    inputSTL = "mycircle.stl" # Enter input STL file name
#    modelName = "solid" # Enter Model Name
#    Writer = MDPAWriter(inputSTL)
#    Writer.CreateModel(modelName)
#    Writer.GiDOutput()





