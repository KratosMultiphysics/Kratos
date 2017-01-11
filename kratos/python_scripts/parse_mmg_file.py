from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *

#gets a line ignoring the empty ones
def GetRelevantLine(f):
    splitted = f.readline().split()
    #print("in reading function " ,splitted)
    while(len(splitted) == 0):
        splitted = f.readline().split()
    return splitted

def ReadVertices(input_file, model_part,parts):
    f = open(input_file,'r')

    #submodel_part_nodes = {}
    #for key,value in parts.items():
        #submodel_part_nodes[key] = []
        
    for line in f:
        splitted = line.split()
        
        if("Vertices" in splitted):
            node_id = 1
            a = GetRelevantLine(f)
            nvertices = int(a[0])
            for i in range(nvertices):
                tmp = GetRelevantLine(f)
                
                #create a model part in the main modelpart
                n = model_part.CreateNewNode(node_id, float(tmp[0]),float(tmp[1]),float(tmp[2]))
                
                #color = int(tmp[3])
                ##print("color = ",color)
                #submodel_part_nodes[color].append(n)
            
                node_id += 1
                
    f.close()
    
    # NOTE: This is not working, reading the nodes of the elements
    ##now add nodes the the corresponding parts
    #for key,value in parts.items():
        #for part in value:
            #if(part != model_part.Name): #nodes are already added to the main model part
                #mp = model_part.GetSubModelPart(part)
                
                #for node in submodel_part_nodes[key]:
                    #mp.AddNode(node,0)
                    
def AddNodesToSubmodelParts(model_part, parts):
    
    for key,value in parts.items():
        for part in value:
            if(part != model_part.Name):
                sub_model_part = model_part.GetSubModelPart(part)
                node_id_list = []
                for elem in sub_model_part.Elements:
                    for node in elem.GetNodes():
                        node_id_list.append(node.Id)
                        
                aux_node_id_list = []
                for i in node_id_list:
                    if i not in aux_node_id_list:
                        aux_node_id_list.append(i)
                    
                sub_model_part.AddNodes(aux_node_id_list)
                    
def ReadTetrahedra(input_file, model_part, parts):
    f = open(input_file,'r')

    submodel_part_elements = {}
    for key,value in parts.items():
        submodel_part_elements[key] = []

    for line in f:
        splitted = line.split()
        
        if("Tetrahedra" in splitted):
            elem_id = 1
            a = GetRelevantLine(f)
            nvertices = int(a[0])
            for i in range(nvertices):
                tmp = GetRelevantLine(f)
                ids = [int(tmp[0]),int(tmp[1]),int(tmp[2]),int(tmp[3])]
                prop_id = int(tmp[4])
                prop = model_part.AddProperties(Properties(prop_id))
                #print(ids)
                e = model_part.CreateNewElement("Element3D4N", elem_id, ids, model_part.GetProperties()[prop_id])
                
                color = int(tmp[4])
                submodel_part_elements[color].append(e)
                
                elem_id += 1
                
    f.close()          
    
    #now add elements the the corresponding parts
    for key,value in parts.items():
        for part in value:
            if(part != model_part.Name): #nodes are already added to the main model part
                mp = model_part.GetSubModelPart(part)
                
                for elem in submodel_part_elements[key]:
                    mp.AddElement(elem,0)    
                    

def ReadTriangles(input_file, model_part,parts):
    f = open(input_file,'r')
    
    submodel_part_conditions = {}
    for key,value in parts.items():
        submodel_part_conditions[key] = []
    
    for line in f:
        splitted = line.split()
        
        if("Triangles" in splitted):
            elem_id = 1
            a = GetRelevantLine(f)
            nvertices = int(a[0])
            for i in range(nvertices):
                tmp = GetRelevantLine(f)
                ids = [int(tmp[0]),int(tmp[1]),int(tmp[2])]
                prop_id = int(tmp[3])
                prop = model_part.AddProperties(Properties(prop_id))
                
                c = model_part.CreateNewCondition("SurfaceCondition3D3N", elem_id, ids, model_part.GetProperties()[prop_id])
                
                color = int(tmp[3])
                submodel_part_conditions[color].append(c)
                
                elem_id += 1
                
    f.close()        
    
    #now add conditions the the corresponding parts
    for key,value in parts.items():
        for part in value:
            if(part != model_part.Name): #nodes are already added to the main model part
                mp = model_part.GetSubModelPart(part)
                
                for cond in submodel_part_conditions[key]:
                    mp.AddCondition(cond,0)    

def ParseMmgFile(input_file, model_part_name = "MainModelPart"):
    ##read kratos input file
    model_part = ModelPart(model_part_name)
    
    #read parts and construct submodelparts
    import json
    with open('part_colors.json') as data_file:    
        tmp = json.load(data_file)
    parts = {}
    for key,value in tmp.items():
        parts[int(key)] = value
        print("Value = ",value)
        if(key != 0): #value of 0 is the main modelpart
            for name in value:
                if(name != model_part.Name):
                    if(not model_part.HasSubModelPart(name)):
                        model_part.CreateSubModelPart(name)
                else:
                    pass #we do not create again the main model part!!
    
    ReadVertices(input_file, model_part, parts)
    
    ReadTriangles(input_file, model_part, parts)
    
    ReadTetrahedra(input_file, model_part, parts)
    
    AddNodesToSubmodelParts(model_part, parts)
    
    print(model_part)
    
    ## NOTE: Esto es una chapuza, arreglar
    #model_part.AddNodalSolutionStepVariable(DISTANCE) 
    #model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT) 
    
    ## NOTE: En C++ poner pnode->SetSolutionStepVariablesList( this_model_part.NodesBegin()->pGetVariablesList() );
    ## Para copiar lista de variables
    #for node in model_part.Nodes:
        ##node.GetSolutionStepValue(DISTANCE, 0)
        #node.SetSolutionStepValue(DISTANCE, 0, abs(node.X))
    
    #from gid_output_process import GiDOutputProcess
    #gid_output = GiDOutputProcess(model_part,
                                #"gid_output",
                                #Parameters("""
                                    #{
                                        #"result_file_configuration" : {
                                            #"gidpost_flags": {
                                                #"GiDPostMode": "GiD_PostAscii",
                                                #"WriteDeformedMeshFlag": "WriteUndeformed",
                                                #"WriteConditionsFlag": "WriteConditions",
                                                #"MultiFileFlag": "SingleFile"
                                            #},        
                                            #"nodal_results"       : []
                                        #}
                                    #}
                                    #""")
                                #)
    #gid_output.ExecuteInitialize()
    #gid_output.ExecuteBeforeSolutionLoop()

    #gid_output.ExecuteInitializeSolutionStep()
    #gid_output.PrintOutput()
    #gid_output.ExecuteFinalizeSolutionStep()
    #gid_output.ExecuteFinalize()   
    
    return model_part
    
if __name__ == '__main__':
    import sys
    model_part = ParseMmgFile(sys.argv[1]) 