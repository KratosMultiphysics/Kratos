from KratosMultiphysics import *

def ReadModelPart(model_part, inputfile):
    #read the json
    import json    
    with open(inputfile,'r') as data_file:    
        data = json.load(data_file)
    
    #read nodes 
    json_nodes = data["model_part"]["Nodes"]
    for item in json_nodes:
        model_part.CreateNewNode(item[0], item[1], item[2], item[3])
    json_nodes = 0 ##free memory not needed anymore
    
    #read properties
    json_properties = data["model_part"]["Properties"]
    for key,val in json_properties.items():
        prop_id = int(key)
        prop = model_part.Properties[prop_id]
        
        variables = val["Variables"]
        for k,v in variables.items():
            var = globals()[str(k)]
            prop[var] = float(v)
    json_properties = 0 ##free memory not needed anymore
    
    #read elements and conditions
    json_elements = data["model_part"]["Elements"]
    for element_name, eldata in json_elements.items():
        connectivity = eldata["connectivity"]
        for tmp in connectivity:
            elem_id = tmp[0]
            prop_id = tmp[1]
            prop = model_part.Properties[prop_id]
            node_ids = []
            for i in range(2,len(tmp)):
                node_ids.append(int(tmp[i]))
            model_part.CreateNewElement(element_name, elem_id, node_ids, prop )
    json_elements = 0 ##free memory not needed anymore
    

    json_conditions = data["model_part"]["Conditions"]
    for condition_name, eldata in json_conditions.items():
        connectivity = eldata["connectivity"]
        for tmp in connectivity:
            elem_id = tmp[0]
            prop_id = tmp[1]
            prop = model_part.Properties[prop_id]
            node_ids = []
            for i in range(2,len(tmp)):
                node_ids.append(int(tmp[i]))
            model_part.CreateNewCondition(condition_name, elem_id, node_ids, prop )
    json_conditions = 0 ##free memory not needed anymore
    
    #read meshes
    json_meshes = data["model_part"]["Meshes"]
    for mesh_name, mesh_data in json_meshes.items():
        #read node pointers
        mesh_id = int(mesh_name)
        mesh = model_part.GetMesh(mesh_id)

        for item in mesh_data["NodePointers"]:
            pnode = model_part.Nodes[item]
            mesh.Nodes.append(pnode)
        
        for item in mesh_data["ElementPointers"]:
            pelem = model_part.Elements[item]
            mesh.Elements.append(pelem)

        for item in mesh_data["ConditionPointers"]:
            pcond = model_part.Conditions[item]
            mesh.Conditions.append(pcond)        
    json_meshes = 0
    
    #read nodal values
    json_nodalvalues = data["model_part"]["NodalData"]
    for varname, vardata in json_nodalvalues.items():
        var = globals()[varname]
        stepindex = vardata["stepindex"]
        vardata = vardata["data"]
        
        for item in vardata:
            node = model_part.Nodes[item[0]]
            if(bool(item[1])):
                node.Fix(var)
            node.SetSolutionStepValue(var,stepindex,float(item[2]))
        
    
    #read conditional values
    
    
    print(model_part)
    