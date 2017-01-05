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

    submodel_part_nodes = {}
    for key,value in parts.items():
        submodel_part_nodes[key] = []
        
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
                
                color = int(tmp[3])
                #print("color = ",color)
                submodel_part_nodes[color].append(n)
            
                node_id += 1
                
    f.close()
    
    
    #now add nodes the the corresponding parts
    for key,value in parts.items():
        for part in value:
            if(part != model_part.Name): #nodes are already added to the main model part
                mp = model_part.GetSubModelPart(part)
                
                for node in submodel_part_nodes[key]:
                    mp.AddNode(node,0)

def ReadTetrahedra(input_file, model_part):
    f = open(input_file,'r')
    
    

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
                model_part.CreateNewElement("Element3D4N", elem_id, ids, model_part.GetProperties()[prop_id])
                elem_id += 1
                
    f.close()            

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
                
                c = model_part.CreateNewCondition("Condition3D3N", elem_id, ids, model_part.GetProperties()[prop_id])
                
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

def ParseMmgFile(input_file):
    ##read kratos input file
    model_part = ModelPart("Main")
    
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
                
    print(model_part)

    
    ReadVertices(input_file, model_part, parts)
    
    ReadTetrahedra(input_file, model_part)
    
    ReadTriangles(input_file, model_part,parts)
    
    print(model_part)
    
    return model_part
    
    
    
if __name__ == '__main__':
    import sys
    model_part = ParseMmgFile(sys.argv[1])

    from gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(model_part,
                                "gid_output",
                                Parameters("""
                                    {
                                        "result_file_configuration" : {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostBinary",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "WriteConditionsFlag": "WriteConditions",
                                                "MultiFileFlag": "SingleFile"
                                            },        
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