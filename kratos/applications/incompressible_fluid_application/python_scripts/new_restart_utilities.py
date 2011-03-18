##def PrintRestart_VectorVariable(variable,varname_x,varname_y,varname_z,Nodes,outfile):
##    outfile.write( "\n" )
##    for node in Nodes:
##        vel = node.GetSolutionStepValue(variable)
##        outstring = "NODES[" + str(node.Id) + "]("+ varname_x + ",0) = "+ str(vel[0]) + ";\n"
##        outstring += "NODES[" + str(node.Id) + "]("+ varname_y + ",0) = "+ str(vel[1]) + ";\n"
##        outstring += "NODES[" + str(node.Id) + "]("+ varname_z + ",0) = "+ str(vel[2]) + ";\n"
##        outfile.write( outstring )

def PrintProperties(outfile):
    outfile.write( "\n" )
    outfile.write( "Begin Properties 0 \n" )
    outfile.write( "End Properties \n" )
    
def PrintNodes(Nodes,outfile):
    outfile.write( "\n" )
    outfile.write( "Begin Nodes \n" )
    for node in Nodes:
        outstring = str(node.Id) + " " + str(node.X)  + " " + str(node.Y) + " " + str(node.Z) + "\n"
        outfile.write( outstring )
    outfile.write( "End Nodes\n" )
    
def PrintElements(element_type,Elements,outfile):
    outfile.write( "\n" )
    outstring = "Begin Elements " + element_type + "\n"
    outfile.write( outstring )
    for elem in Elements:
        outstring = str(elem.Id) + " " + "0" # str(elem.Properties.Id())
        for node in elem.GetNodes():
            outstring += " " + str(node.Id)
        outstring += "\n"
        outfile.write( outstring )
    outfile.write( "End Elements\n" )
    
def PrintConditions(element_type,Conditions,outfile):
    outfile.write( "\n" )
    outstring = "Begin Conditions " + element_type + "\n"
    outfile.write( outstring )
    for elem in Conditions:
        outstring = str(elem.Id) + " " + "0" # str(elem.Properties.Id())
        for node in elem.GetNodes():
            outstring += " " + str(node.Id)
        outstring += "\n"
        outfile.write( outstring )
    outfile.write( "End Conditions\n" )
    
def PrintRestart_ScalarVariable(variable,varname,Nodes,outfile):
    outfile.write( "\n" )
    outstring = "Begin NodalData " + varname + "\n"
    outfile.write( outstring )
    for node in Nodes:
        is_fixed = 0
        if(node.IsFixed(variable) ):
            is_fixed = 1              
        outstring = str(node.Id) + " " + str(is_fixed) + " "+ str(node.GetSolutionStepValue(variable)) + "\n"
        outfile.write( outstring )
    outfile.write( "End NodalData\n" )   



