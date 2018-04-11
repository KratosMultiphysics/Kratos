from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *

white_listed_vars = (
                     "",
                     "FILLTIME (s)", 
                     "LAST AIR", 
                     "MAX VEL (m/s)",
                     "PRESSURES (N/m2)",
                     "SHRINKAGE_POROSITY (in^3)",
                     "SHRINKAGE_POROSITY (m^3)",
                     "SOLID FRACTION",
                     "SOLIDIF MODULUS (cm)",
                     "SOLIDIF MODULUS (in)",
                     "SOLIDIF TIME (s)"
                     )

def PrintVariablesString(header_name,input_string):
    input_string = input_string.replace("    ","")
    ##
    input_string = input_string.rstrip()
    
    #ititialize block
    
    words = input_string.split('\n')
    
    outstring = str("")
    
    for word in words:
        #print("verifying variable ",word, )
        if(globals().get(word) == None): 
            if word not in white_listed_vars:
                outstring += str("WARNING: variable ") + str(word) + str(" is not correctly exported to python -->key is: ") + str( globals().get(word))  + str("\n")
        if(globals().get(word) == 0): 
            if word not in white_listed_vars:
                outstring += str("WARNING: variable ") + str(word) + str(" is not correctly exported to python -->key is: ") + str( globals().get(word))  + str("\n")

    return outstring
        
    
def CheckVariables():
    
    kernel = KratosGlobals.Kernel
    
    
    #print double Variables
    var_list = kernel.GetDoubleVariableNames()
    var_list += kernel.GetVariableComponentVariableNames()
    
    warning_list = PrintVariablesString("KratosVar_double",var_list)
        
    #print integer variables
    var_list = kernel.GetIntVariableNames()
    var_list += kernel.GetUnsignedIntVariableNames() 
    warning_list += PrintVariablesString("KratosVar_int",var_list)
        
    #print Bool Variables
    var_list = kernel.GetBoolVariableNames()
    warning_list += PrintVariablesString("KratosVar_bool",var_list)
    
    #print array_1d Variables
    var_list = kernel.GetArrayVariableNames()
    warning_list += PrintVariablesString("KratosVar_Array1d",var_list)

    #print Vector Variables
    var_list = kernel.GetVectorVariableNames()
    warning_list += PrintVariablesString("KratosVar_Vector",var_list)
    
    #print Matrix Variables
    var_list = kernel.GetMatrixVariableNames()
    warning_list += PrintVariablesString("KratosVar_Matrix",var_list)
    
    #print Flags
    var_list = kernel.GetFlagsVariableNames()
    warning_list += PrintVariablesString("KratosVar_FlagsType",var_list)
    
    return warning_list

print(CheckVariables())