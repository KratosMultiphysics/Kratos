from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.MeshingApplication import *



def PrintVariablesString(header_name,input_string):
    #ititialize block
    outstring = "\n<xs:simpleType name=\"" + header_name + "\"> \n"
    outstring += "    <xs:restriction base=\"xs:string\"> \n"
    
    words = input_string.split()
    
    for word in words:
        outstring += "         <xs:enumeration value=\"" + word + "\"> \n"
        
    #finalize block
    outstring += "     </xs:restriction> \n"
    outstring += "</xs:simpleType> \n"
    outstring += "\n"
    
    return outstring
    
def PrintVariablesToXMLenum(filename):
    
    kernel = KratosGlobals.Kernel
    
    out = open(filename,'w')
    
    #print XSD header
    out.write(" <?xml version=\"1.0\"?> \n")
    out.write("<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\" elementFormDefault=\"qualified\"> \n")
    
    #print double Variables
    var_list = kernel.GetDoubleVariableNames()
    tmp = PrintVariablesString("KratosVar_double",var_list)
    out.write(tmp)
        
    #print integer variables
    var_list = kernel.GetIntVariableNames()
    tmp = PrintVariablesString("KratosVar_int",var_list)
    out.write(tmp)
    
    #print Bool Variables
    var_list = kernel.GetBoolVariableNames()
    tmp = PrintVariablesString("KratosVar_bool",var_list)
    out.write(tmp)    
    
    #print array_1d Variables
    var_list = kernel.GetArrayVariableNames()
    tmp = PrintVariablesString("KratosVar_Array1d",var_list)
    out.write(tmp)      
    
    #finalize XSD schema and close file
    out.write("</xs:schema")
    out.close()


PrintVariablesToXMLenum("vars.xmd")