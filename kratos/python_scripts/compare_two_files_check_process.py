from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

import filecmp 
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CompareTwoFilesCheckProcess(Model, settings["Parameters"])

class CompareTwoFilesCheckProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
  
    def __init__(self,model_part,params):
        
        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "file_name_1"            : "",
            "file_name_2"            : "",
            "deterministic"          : true,
            "non_deterministic_comp" : "mesh_file",
            "error_assumed"          : 1.0e-6,
            "dimension"              : 3
        }
        """)
        
        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
        self.dimension = self.params["dimension"].GetInt()
        self.non_deterministic_comp = self.params["non_deterministic_comp"].GetString()
        self.deterministic = self.params["deterministic"].GetBool()
        self.error_assumed = self.params["error_assumed"].GetDouble()
        
    def ExecuteInitialize(self):
        self.file_name_1 = os.path.join(os.getcwd(), self.params["file_name_1"].GetString())
        self.file_name_2 = os.path.join(os.getcwd(), self.params["file_name_2"].GetString())
        
    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass
            
    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        if (self.deterministic == True):
            value = filecmp.cmp(self.file_name_1, self.file_name_2)
            self.assertTrue(value)
        else:
            if (self.non_deterministic_comp == "mesh_file"):
                error = _ReadVertices(self.file_name_1, self.file_name_2, self.dimension)
                self.assertTrue(error < self.error_assumed)
            elif (self.non_deterministic_comp == "sol_file"):
                error = _ReadMetric(self.file_name_1, self.file_name_2, self.dimension)
                self.assertTrue(error < self.error_assumed)
            else: 
                raise NameError('Non-deterministic comparision not implemented yet')

def _ConvertStringToListFloat(line, space = " ", endline = ""):
    list_values = []
    string_values = (line.replace(endline,"")).split(space)
    for string in string_values:
        list_values.append(float(string))
    #value = ""
    #for i in line:
        #if ((i == " ") or (i == "  ")):
            #num_value = float(value)
            #list_values.append(num_value)
            #value = ""
        #else:
            #value += i
    return list_values

def _ReadVertices(input_file1, input_file2, dimension):
    f1 = open(input_file1,'r')
    f2 = open(input_file2,'r')
    
    lines1 = f1.readlines()
    lines2 = f2.readlines()
    
    numline = 0
    for line1 in lines1:
        numline += 1
        
        if("Vertices" in line1):
            line = lines1[numline]
            nvertices = int(line)
            numline += 1
            break
        
    error = 0.0
    for i in range(numline, nvertices + numline):
        tmp1 = _ConvertStringToListFloat(lines1[i], "", "\n")
        tmp2 = _ConvertStringToListFloat(lines2[i], "", "\n")
        if (dimension == 2): 
            error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0)**(0.5)
        else:
            error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0 + (tmp1[2] - tmp2[2])**2.0)**(0.5)    
            
    f1.close()
    f2.close()
    
    return (error/nvertices)

def _ReadMetric(input_file1, input_file2, dimension):
    f1 = open(input_file1,'r')
    f2 = open(input_file2,'r')
    
    lines1 = f1.readlines()
    lines2 = f2.readlines()
    
    numline = 0
    for line1 in lines1:
        numline += 1
        
        if("SolAtVertices" in line1):
            line = lines1[numline]
            nvertices = int(line)
            numline += 2
            break
        
    error = 0.0
    for i in range(numline, nvertices + numline):
        
        if dimension == 2:
            space = " "
            end_line = " \n"
        else:
            space = "  "
            end_line = "  \n"
        
        if (lines1[i][0] == " "):
            lines1[i] = lines1[i][1:]
        if (lines2[i][0] == " "):
            lines1[i][0] = lines2[i][1:]
        tmp1 = _ConvertStringToListFloat(lines1[i], space, end_line)
        tmp2 = _ConvertStringToListFloat(lines2[i], space, end_line)
        
        if (dimension == 2): 
            error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0 + (tmp1[2] - tmp2[2])**2.0)**(0.5)   
        else:
            error += ((tmp1[0] - tmp2[0])**2.0 + (tmp1[1] - tmp2[1])**2.0 + (tmp1[2] - tmp2[2])**2.0 + (tmp1[3] - tmp2[3])**2.0 + (tmp1[4] - tmp2[4])**2.0 + (tmp1[5] - tmp2[5])**2.0)**(0.5)    
            
    f1.close()
    f2.close()
    
    return (error/nvertices)
