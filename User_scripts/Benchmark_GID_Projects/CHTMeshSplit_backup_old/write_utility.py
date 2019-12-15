from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

def PrintModelData(outfile):
    outfile.write("\n")
    outfile.write("Begin ModelPartData \n")
    outfile.write("End ModelPartData \n")

def PrintProperties(outfile):
    outfile.write("\n")
    outfile.write("Begin Properties 0 \n")
    outfile.write("End Properties \n")


def PrintNodes(Nodes, outfile):
    outfile.write("\n")
    outfile.write("Begin Nodes \n")
    for node in Nodes:
        outstring = str(node.Id) + " " + str(node.X) + " " + str(node.Y) + " " + str(node.Z) + "\n"
        outfile.write(outstring)
    outfile.write("End Nodes\n")


def PrintElements(element_type, Elements, outfile):
    outfile.write("\n")
    outstring = "Begin Elements " + element_type + "\n"
    outfile.write(outstring)
    for elem in Elements:
        outstring = str(elem.Id) + "   " + str(elem.Properties.Id)
        for node in elem.GetNodes():
            outstring += "   " + str(node.Id)
        outstring += "\n"
        outfile.write(outstring)
    outfile.write("End Elements\n")


def PrintConditions(element_type, Conditions, outfile):
    outfile.write("\n")
    outstring = "Begin Conditions " + element_type + "\n"
    outfile.write(outstring)
    for elem in Conditions:
        outstring = str(elem.Id) + "   " + str(elem.Properties.Id)
        for node in elem.GetNodes():
            outstring += "   " + str(node.Id)
        outstring += "\n"
        outfile.write(outstring)
    outfile.write("End Conditions\n")


def PrintRestart_ScalarVariable(variable, varname, Nodes, outfile):
    outfile.write("\n")
    outstring = "Begin NodalData " + varname + "\n"
    outfile.write(outstring)
    for node in Nodes:
        is_fixed = 0
        if(node.IsFixed(variable)):
            is_fixed = 1
        outstring = str(
            node.Id) + " " + str(
                is_fixed) + " " + str(
                    node.GetSolutionStepValue(
                        variable)) + "\n"
        outfile.write(outstring)
    outfile.write("End NodalData\n")

def PrintRestart_Variable_On_Condition(variable, varname, Conditions, outfile):
    outfile.write("\n")
    outstring = "Begin ConditionalData " + varname + "\n"
    outfile.write(outstring)
    for cond in Conditions:
        outstring = str(cond.Id) + " "  + str(cond.GetValue(variable)) + "\n"
        outfile.write(outstring)
    outfile.write("End ConditionalData\n")