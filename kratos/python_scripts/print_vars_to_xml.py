from KratosMultiphysics import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *



def PrintVariablesString(header_name,input_string):
    input_string = input_string.replace("    ","")
    #ititialize block
    outstring = "\n<xs:simpleType name=\"" + header_name + "\"> \n"
    outstring += "    <xs:restriction base=\"xs:string\"> \n"

    words = input_string.split('\n')

    for word in words:
        outstring += "         <xs:enumeration value=\"" + word + "\"> \n"

    #finalize block
    outstring += "     </xs:restriction> \n"
    outstring += "</xs:simpleType> \n"
    outstring += "\n    "

    return outstring

def PrintVariablesToXMLenum(filename):

    kernel = KratosGlobals.Kernel

    out = open(filename,'w')

    #print XSD header
    out.write(" <?xml version=\"1.0\"?> \n")
    out.write("<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\" elementFormDefault=\"qualified\"> \n")

    #print double Variables
    var_list = kernel.GetDoubleVariableNames()
    var_list += kernel.GetVariableComponentVariableNames()

    tmp = PrintVariablesString("KratosVar_double",var_list)
    out.write(tmp)

    #print integer variables
    var_list = kernel.GetIntVariableNames()
    var_list += kernel.GetUnsignedIntVariableNames()
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

    #print Vector Variables
    var_list = kernel.GetVectorVariableNames()
    tmp = PrintVariablesString("KratosVar_Vector",var_list)
    out.write(tmp)

    #print Matrix Variables
    var_list = kernel.GetMatrixVariableNames()
    tmp = PrintVariablesString("KratosVar_Matrix",var_list)
    out.write(tmp)

    #print Flags
    var_list = kernel.GetFlagsVariableNames()
    tmp = PrintVariablesString("KratosVar_FlagsType",var_list)
    out.write(tmp)

    #finalize XSD schema and close file
    out.write("</xs:schema")
    out.close()


PrintVariablesToXMLenum("kratos_variables.xsd")