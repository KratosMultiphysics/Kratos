import os, sys
from time import *

def WriteInfo(i_file, write_mode, info):
    file_for_output = open(i_file ,write_mode)
    file_for_output.write("=== " + info + " ===\n")
    file_for_output.flush()
    file_for_output.close()
    print(info)

input_file = "MainKratos_Mapping_Parallel.py"
output_file = "kratos_output.txt"

list_processors = []

for i in range(1,31):
    list_processors.append(i)

list_processors.extend([40, 50])

# serial execution
WriteInfo(output_file, "w", "Serial Execution")
system_cmd = "python3.5 " + input_file + " >> " + output_file
os.system(system_cmd)

# parallel executions
for num_processors in list_processors:
    WriteInfo(output_file, "a", "Parallel Execution; " + str(num_processors) + " processors")
    system_cmd = "mpiexec -np " + str(num_processors) + " python3.5 " + input_file + " >> " + output_file
    os.system(system_cmd)
