# Philipp Bucher, 2.12.2016
# This script is intended to test the MappingApplication in parallel with
# several numbers of cores, since at the time of this work there is no parallel
# testing in Kratos yet
# usage: "python3 runscript_test_MappingApplication.py"

import os, sys
from time import *
import datetime
import re

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

def CheckOutputFile(output_file, string_array, tests_success):
    with open(output_file, 'r') as test_file:
        i = 0
        for line in test_file:
            i += 1
            for keyword in string_array:
                if re.search(keyword, line, re.IGNORECASE):
                    tests_success = False
                    print(keyword + " in line " + str(i) + " of file " + output_file)

    return tests_success

def DeleteAuxTestFiles(Path, MDPAFileNames):
    # This function deletes the time and the partitioned files
    files = os.listdir(Path)
    for file in files:
        if file not in MDPAFileNames:
            os.remove(os.path.join(Path,file))


input_file = "test_MappingApplication.py"
kratos_output_file = "output_kratos.txt"
tests_output_file = "output_tests.txt"

list_processors = []

for i in range(1,31):
    list_processors.append(i)

list_processors.extend([40, 50])

os.system("export OMP_NUM_THREADS=1")
print("OMP thread set to 1")

# serial execution
WriteInfo(kratos_output_file, tests_output_file, "w", "Serial Execution")
system_cmd = "python3 " + input_file + " >> " + kratos_output_file + " 2>> " + tests_output_file
os.system(system_cmd)

# parallel executions
for num_processors in list_processors:
    WriteInfo(kratos_output_file, tests_output_file, "a", "Parallel Execution; " + str(num_processors) + " processors")
    system_cmd = "mpiexec -np " + str(num_processors) + " python3 " + input_file_mpi + " >> " + kratos_output_file + " 2>> " + tests_output_file
    os.system(system_cmd)

tests_success = True

keyword_array = ["FAIL", "mpiexec", "mpirun", "Segmentation", "signal", "not"]
keyword_array.extend(["Traceback", "RuntimeError", "ERROR", "Error", "WARNING", "Errno"])

tests_success = CheckOutputFile(tests_output_file, keyword_array, tests_success)
tests_success = CheckOutputFile(kratos_output_file, keyword_array, tests_success)

test_runtime = datetime.timedelta(seconds=int((time() - start_time)))
print("\nTests Sucessful: " + str(tests_success) + ", Runtime: " + str(test_runtime))


## Removing the files that are created during the tests (patritioned mdpas and time files)
tests_execution_path = os.path.dirname(os.path.realpath(__file__))
mdpa_files_path = os.path.join(tests_execution_path, "MapperTests_mdpa")
DeleteAuxTestFiles(mdpa_files_path,
                   ["MappingApplication_test_geometry_quad.mdpa",
                    "MappingApplication_test_geometry_tri.mdpa"])

mdpa_files_path = os.path.join(tests_execution_path, "NearestNeighborMapperTest_mdpa")
DeleteAuxTestFiles(mdpa_files_path,
                   ["origin.mdpa",
                    "destination.mdpa"])

mdpa_files_path = os.path.join(tests_execution_path, "NearestElementMapperTest2D_mdpa")
DeleteAuxTestFiles(mdpa_files_path,
                   ["origin.mdpa",
                    "destination.mdpa"])
