from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *

def ComputeTensorH(d,v,ratio,h):
    if(d>h):
        coeff = 1/(h*h)
        mystring = str(coeff)+" 0 "+str(coeff) + " 0 0 "+str(coeff) + str("\n")
    else:
        
        coeff1 = 1/h*h
        
        a = ratio + (d/h)*(1-ratio)

        coeff2 = coeff1/(a*a)
        mystring = str(coeff1*(1-v[0]*v[0]) + coeff2*v[0]*v[0]) + " "
        mystring += str(coeff1*(v[0]*v[1]) + coeff2*v[0]*v[1]) + " "
        mystring += str(coeff1*(1-v[1]*v[1]) + coeff2*v[1]*v[1]) + " "
        mystring += str(coeff1*(v[0]*v[2]) + coeff2*v[0]*v[2]) + " "
        mystring += str(coeff1*(v[1]*v[2]) + coeff2*v[1]*v[2]) + " "
        mystring += str(coeff1*(1-v[2]*v[2]) + coeff2*v[2]*v[2]) + " "
        mystring +=" \n"
    return mystring


##here we assign a color to all the nodes and conditions so to assign them a "reference" and be able to reassign them to the correct part
##TODO: rough around the edges. only works with one level of subparts
def ComputeColors(input_file,model_part):
    
    
    ###TODO: remove just for testing
    #sub1 = model_part.CreateSubModelPart("sub1")
    #sub2 = model_part.CreateSubModelPart("sub2")    
    #sub1.AddNode(model_part.Nodes[1],0)
    #sub2.AddNode(model_part.Nodes[2],0)

    
    submodel_colors = {}
    
    #obtain a flat list of submodelparts:
    parts = []
    parts.append(model_part)
    for subpart in model_part.SubModelParts:
        name = subpart.Name 
        parts.append(model_part.GetSubModelPart(name))
        for subsub in subpart.SubModelParts:
            name = subpart.Name
            parts.append(subsub.GetSubModelPart(name))
        
    node_colors = {}
    cond_colors = {}
    for node in model_part.Nodes:
        node_colors[node.Id] = set()
    for cond in model_part.Conditions:
        cond_colors[cond.Id] = set()
        
    color = 0
    for part in parts:
        submodel_colors[color] = [part.Name]
        
        if(color != 0): #general model part would overlap with everythign
            for node in part.Nodes:
                node_colors[node.Id].add(color)
            for cond in part.Conditions:
                cond_colors[cond.Id].add(color)
                
        color+=1
        
        
    #now detect all the cases in which a node or a cond belongs to more than one part simultaneously
    combinations = {}
    for key,value in node_colors.items():
        if(len(value) > 1):
            combinations[frozenset(value)] = -1
    for key,value in cond_colors.items():
        if(len(value) > 1):
            combinations[frozenset(value)] = -1   
            
    for key,value in combinations.items():
        submodel_colors[color] = []
        #print("key =",key,value)
        for item in key:
            #print("item = ", item)
            submodel_colors[color].append(  submodel_colors[item][0] )
        combinations[key] = color
        color += 1
        
    print("combinations = ",combinations)
        
    #write equivalence to a file
    import json
    with open('part_colors.json', 'w') as f:
        json.dump(submodel_colors, f, ensure_ascii=False)
        
    #now overwrite the node_colors and cond_colors with the corresponding value
    for key,value in node_colors.items():
        tmp = value
        if(len(value) == 0):
            node_colors[key] = 0 #it belongs to the main model part
        elif(len(value)==1): #get the value
            node_colors[key] = list(tmp)[0]
        else:
            node_colors[key] = combinations[frozenset(tmp)]

    for key,value in cond_colors.items():
        tmp = value
        if(len(value) == 0):
            cond_colors[key] = 0 #it belongs to the main model part
        elif(len(value)==1): #get the value
            cond_colors[key] = list(tmp)[0]
        else:
            print()
            cond_colors[key] = combinations[frozenset(tmp)]
    #print("node_colors = ",node_colors)
    return node_colors,cond_colors

        
    


def WriteMmgFile(input_file):
    ##read kratos input file
    model_part = ModelPart("Main")
            
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)

    ModelPartIO(input_file).ReadModelPart(model_part)
    
    ##ensure ordering is consecutive
    index = 1
    for node in model_part.Nodes:
        node.Id = index
        index+=1
    index = 1
    for cond in model_part.Conditions:
        cond.Id = index
        index+=1
    index = 1
    for elem in model_part.Elements:
        elem.Id = index
        index+=1
        
    node_colors,cond_colors= ComputeColors(input_file, model_part)
    
    
    ##open mmg output file
    output_file = input_file+".mesh"
    outfile = open(output_file,'w')

    outfile.write("MeshVersionFormatted 2\n\n")
    outfile.write("Dimension 3\n\n")

    #write nodes
    outfile.write("Vertices\n")
    outfile.write(str(len(model_part.Nodes))+"\n")
    for node in model_part.Nodes:
        outline = str(node.X)+" "+str(node.Y) + " " +str(node.Z)+ " " + str(node_colors[node.Id]) +"\n"
        outfile.write(outline)
    outfile.write("\n")

    #write Tetras
    if(len(model_part.Elements) != 0):
        outfile.write("Tetrahedra\n")
        outfile.write(str(len(model_part.Elements))+"\n")
        for elem in model_part.Elements:
            outline = ""
            for node in elem.GetNodes():
                outline += str(node.Id) + " " 
            outline += "1\n"
            outfile.write(outline)
        outfile.write("\n")
        outfile.write("\n")
    
    #write Triangles
    if(len(model_part.Conditions) != 0):
        outfile.write("Triangles\n")
        outfile.write(str(len(model_part.Conditions))+"\n")
        for cond in model_part.Conditions:
            outline = ""
            for node in cond.GetNodes():
                outline += str(node.Id) + " " 
            outline += str(cond_colors[cond.Id]) + "\n"
            outfile.write(outline)
        outfile.write("\n")

    outfile.write("End\n")
    outfile.close()
    
    #print("bbb")
    ####here write the .sol file
    sol_file = input_file+".sol"
    sol_file = open(sol_file,'w')

    sol_file.write("MeshVersionFormatted 2\n\n")
    sol_file.write("Dimension 3\n\n")
    sol_file.write("SolAtVertices\n")
    sol_file.write(str(len(model_part.Nodes))+"\n")

    sol_file.write("1 1\n") #scalar
    for i in range(len(model_part.Nodes)):
        sol_file.write(str(0.0010)+"\n")
        #sol_file.write("200.0 0 200.0 0 0 10000.0\n")
        
    #sol_file.write("1 3\n") #tensor
    #v = Vector(3)
    #v[0] = 1.0
    #v[1] = 0.0
    #v[2] = 0.0
    #h = 0.1
    #ratio = 0.01
    #for node in model_part.Nodes:
        #d = node.X*0.1
        #sol_file.write(ComputeTensorH(d,v,ratio,h))
        
        
    sol_file.write("End\n")
    sol_file.close()
    
    print("finished")
    
if __name__ == '__main__':
    import sys
    WriteMmgFile(sys.argv[1])