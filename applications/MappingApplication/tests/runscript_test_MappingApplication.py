# Philipp Bucher, 2.12.2016
# This script is intended to test the Nearest Neighbor Mapper in parallel with 
# several numbers of cores, since at the time of this work there is no parallel
# testing in Kratos yet
# usage: "python3.5 runscript_test_MappingApplication.py"

import os, sys
from time import *
import datetime

start_time = time()

def WriteInfo(kratos_file, test_file, write_mode, info):
    file_kratos_output = open(kratos_file ,write_mode)
    file_kratos_output.write("=== " + info + " ===\n")
    file_kratos_output.flush()
    file_kratos_output.close()

    file_test_output = open(test_file ,write_mode)
    file_test_output.write("\n\n\n=== " + info + " ===\n")
    file_test_output.flush()
    file_test_output.close()

    print(info)

input_file = "test_MappingApplication.py"
kratos_output_file = "output_kratos.txt"
tests_output_file = "output_tests.txt"

list_processors = []

for i in range(1,31):
    list_processors.append(i)

list_processors.extend([40, 50])

# serial execution
WriteInfo(kratos_output_file, tests_output_file, "w", "Serial Execution")
system_cmd = "python3.5 " + input_file + " >> " + kratos_output_file + " 2>> " + tests_output_file
os.system(system_cmd)

# parallel executions
for num_processors in list_processors:
    WriteInfo(kratos_output_file, tests_output_file, "a", "Parallel Execution; " + str(num_processors) + " processors")
    system_cmd = "mpiexec -np " + str(num_processors) + " python3.5 " + input_file + " >> " + kratos_output_file + " 2>> " + tests_output_file
    os.system(system_cmd)

tests_success = True
with open(tests_output_file, 'r') as test_file:
    i = 0
    for line in test_file:
        i += 1
        if "FAIL" in line:
            tests_success = False
            print("FAIL in line " + str(i))

test_runtime = datetime.timedelta(seconds=int((time() - start_time)))
print("\nTests Sucessful: " + str(tests_success) + ", Runtime: " + str(test_runtime))
