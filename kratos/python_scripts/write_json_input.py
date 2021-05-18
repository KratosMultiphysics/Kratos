from KratosMultiphysics import *
import json

class JsonWriter:
    def __init__(self, filename, model_part):
        self.filename = filename
        self.model_part = model_part
        self.json_root = {}
        self.json_root["model_part"] = {}
        
    def DumpConnectivity(self, element_name, condition_name):
        
        model_part_root = self.json_root["model_part"]
    
        model_part_root["PartName"] = self.model_part.Name 
        
        model_part_root["Nodes"] = []
        json_nodes = model_part_root["Nodes"]
        for node in self.model_part.Nodes:
            json_nodes.append([ int(node.Id), float(node.X), float(node.Y), float(node.Z) ] ) 
        
        model_part_root["Elements"] = {element_name : {"connectivity" : [] } }
        connectivity = model_part_root["Elements"][element_name]["connectivity"]
        for elem in self.model_part.Elements:
            tmp = [ elem.Id, elem.Properties.Id ]
            for node in elem.GetNodes():
                tmp.append( node.Id )
            connectivity.append(tmp)
            
        #add elements to the json - assume they are of type "Stokes3D"
        model_part_root["Conditions"] = {condition_name: {"connectivity" : [] } }
        connectivity = model_part_root["Conditions"][condition_name]["connectivity"]
        for elem in self.model_part.Conditions:
            tmp = [ elem.Id, elem.Properties.Id ]
            for node in elem.GetNodes():
                tmp.append( node.Id )
            connectivity.append(tmp)

    def DumpProperties(self, property_variable_list): #, property_table_list):
        model_part_root = self.json_root["model_part"]
        model_part_root["Properties"] = {}
        properties_json_root = model_part_root["Properties"]
        for prop in self.model_part.Properties:
            properties_json_root[ prop.Id ] = {"Variables": {}, "Tables": {} }
            
            for varname in property_variable_list:
                var = globals().get(varname)
                properties_json_root[ prop.Id ]["Variables"][varname] = prop[var]
                
            #for table_name in property_table_list:
                #properties_json_root[ prop.Id ]["Tables"][table_name] = prop.GetTable(table_name)
                
            
    def DumpSolutionStepNodalData(self, nodal_variables_list, stepindex):
        model_part_root = self.json_root["model_part"]
        model_part_root["NodalData"] = {}
        for varname in nodal_variables_list:
            var = globals().get(varname)
            model_part_root["NodalData"][varname] = {"stepindex":stepindex, "data" : []}
            aux = model_part_root["NodalData"][varname]["data"]
            for node in self.model_part.Nodes:
                aux.append([node.Id, int(node.IsFixed(var)), float(node.GetSolutionStepValue(var, stepindex))])
                            
    def DumpMeshes(self):
        model_part_root = self.json_root["model_part"]
        model_part_root["Meshes"] = {}
        
        counter = 0
        for mesh in self.model_part.GetListOfMeshes():
            counter += 1
            meshname = str(counter) #TODO: change this to get the name of the mesh when available
            
            model_part_root["Meshes"][meshname] = {"NodePointers":[], "ElementPointers":[], "ConditionPointers":[]}
            
            #mesh = self.model_part.GetMesh(meshname)
               
            nodepointers = model_part_root["Meshes"][meshname]["NodePointers"]
            for node in mesh.Nodes:
                nodepointers.append( node.Id )
                
            elementpointers = model_part_root["Meshes"][meshname]["ElementPointers"]
            for elem in mesh.Elements:
                elementpointers.append( elem.Id )
                
            conditionpointers = model_part_root["Meshes"][meshname]["ConditionPointers"]
            for cond in mesh.Conditions:
                conditionpointers.append( cond.Id )
                
    def WriteJson(self):
         outfile = open(self.filename, 'w')
         
         print(self.filename)
         
         import json
         json.dump(self.json_root, outfile, indent=4)
         
         outfile.close()
