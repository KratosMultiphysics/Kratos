##def PrintRestart_VectorVariable(variable,varname_x,varname_y,varname_z,Nodes,outfile):
##    outfile.write( "\n" )
##    for node in Nodes:
##        vel = node.GetSolutionStepValue(variable)
##        outstring = "NODES[" + str(node.Id) + "]("+ varname_x + ",0) = "+ str(vel[0]) + ";\n"
##        outstring += "NODES[" + str(node.Id) + "]("+ varname_y + ",0) = "+ str(vel[1]) + ";\n"
##        outstring += "NODES[" + str(node.Id) + "]("+ varname_z + ",0) = "+ str(vel[2]) + ";\n"
##        outfile.write( outstring )

def PrintRestart_ScalarVariable(variable,varname,Nodes,outfile):
    outfile.write( "\n" )
    for node in Nodes:
        temp = node.GetSolutionStepValue(variable)
        outstring = "NODES[" + str(node.Id) + "]("+ varname + ",0) = "+ str(temp) + ";\n"
        outfile.write( outstring )



def PrintRestartFixity(variable,varname,Nodes,outfile):
    outfile.write( "\n" )
    for node in Nodes:
        if( node.IsFixed(variable) ):
            outstring = "NODES[" + str(node.Id) + "].Fix("+ varname + ");\n"
            outfile.write( outstring )
        else:
            outstring = "NODES[" + str(node.Id) + "].Free("+ varname + ");\n"
            outfile.write( outstring )

            



def PrintRestart_ScalarVariable_PyFormat(variable,varname,Nodes,outfile):
    outfile.write( "\n" )
    for node in Nodes:
        temp = node.GetSolutionStepValue(variable)
        outstring = "NODES[" + str(node.Id) + "].SetSolutionStepValue("+ varname + ",0," + str(temp) +");\n"
        outfile.write( outstring )

def PrintRestartFixity_PyFormat(variable,varname,Nodes,outfile):
    outfile.write( "\n" )
    for node in Nodes:
        if( node.IsFixed(variable) ):
            outstring = "NODES[" + str(node.Id) + "].Fix("+ varname + ");\n"
            outfile.write( outstring )
        else:
            outstring = "NODES[" + str(node.Id) + "].Free("+ varname + ");\n"
            outfile.write( outstring )
