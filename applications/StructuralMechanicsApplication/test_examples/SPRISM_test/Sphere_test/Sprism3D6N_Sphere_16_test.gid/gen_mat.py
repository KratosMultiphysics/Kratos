import os

file = open("materials.py", "wb")

file.write(("# Importing the Kratos Library\n").encode('utf-8'))
file.write(("from KratosMultiphysics import *\n").encode('utf-8'))
file.write(("from KratosMultiphysics.SolidMechanicsApplication import *\n").encode('utf-8'))
file.write(("def AssignMaterial(Properties):\n").encode('utf-8'))

cwd_full_path = os.getcwd()  
cwd_dir_name = os.path.basename(cwd_full_path)

lines = tuple(open(cwd_dir_name[:-4]+".mdpa", 'r'))

count = 0
for i in range(len(lines)):
    if "Begin Properties" in lines[i]:
        count+=1
        file.write(("    prop_id = "+str(count)+";\n").encode('utf-8'))
        file.write(("    prop = Properties["+str(count)+"]\n").encode('utf-8'))
        const_law = lines[i + 1][23:-1]
        file.write(("    mat = "+const_law+"();\n").encode('utf-8'))
        file.write(("    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());\n").encode('utf-8'))

file.close() 
