import os
from os.path import join
import subprocess
import time
import numpy as np
import copy


import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure

def FORT_replace(line, orig, new):
    '''The length of a line in FORTRAN 77 is limited, replacing working directories can exceed this limiet
    This functions splits these strings over multiple lines'''

    ampersand_location = 6
    char_limit = 72

    if "|" in line:
        temp = line.replace(orig, new)
        N = len(temp)

        if N > char_limit:
            count = 0
            line = ""
            line += temp[0:char_limit]+"\n"
            count +=char_limit
            while count < N:
                print(count)
                temp_string = temp[count:count+char_limit-12]
                n = len(temp_string)
                count +=n
                if count < N: #need to append an additional new line
                    line+= "     &"+"      "+temp_string+"\n"
                else:
                    line += "     &" + "      " + temp_string
                print (line)
        else:
            line = temp

    return line


# prepare Abaqus USRInit.f
usr = "../../USRInit.f"
with open(usr, "r")as infile:
    with open( "usr.f", "w") as outfile:
        for line in infile:
            line = line.replace("|dimension|", str(3))
            line = line.replace("|surfaces|", str(2))
            line = line.replace("|cpus|", str(1))

            # if PWD is too ling then FORTRAN code can not compile so this needs special treatment
            line = FORT_replace(line, "|PWD|", os.path.abspath(os.path.join("Random_Directory", os.pardir)))
            line = FORT_replace(line, "|CSM_dir|", "Random_Directory/")

            if "|" in line:
                raise ValueError(
                    f"The following line in USRInit.f still contains a \"|\" after substitution: \n \t{line} \n Probably a parameter was not subsituted")
            outfile.write(line)