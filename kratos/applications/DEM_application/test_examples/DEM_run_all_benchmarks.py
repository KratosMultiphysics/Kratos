from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys
import platform

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
dem_scripts_path = '../test_examples/basic_benchmarks'
sys.path.append(dem_scripts_path)
import benchmarking
os.chdir(dem_scripts_path)


def Run():
    
    print("\nStarting DEM Benchmarks..............\n")
    f = open("errors.txt", "w")
    f.write("\n========== DEM BENCHMARKS ===========\n")
    f.write("========== SLIDING REGIME ===========\n\n")
    f.close()
    Text=""
        
    for benchmark in range(3, 9):
      
        print("Running Benchmark " + str(benchmark) + " of 8.............")    
        
        if platform.system()=="Windows":
            os.system("python Chung_Ooi_benchmarks.py " + str(benchmark) + " > BenchTemp.txt")
        else:
            if sys.version_info >= (3, 0):
                os.system("python3 Chung_Ooi_benchmarks.py " + str(benchmark) + " > BenchTemp.txt")
            else:
                os.system("python -3 Chung_Ooi_benchmarks.py " + str(benchmark) + " > BenchTemp.txt")
                
    os.remove("BenchTemp.txt")
        
    f = open("errors.txt")
    file_contents = f.read()
    f.close()
    
    Text += file_contents.rstrip("\n")
    Text += "\n\n\n"
        
    os.remove("errors.txt")
    
    return Text


if __name__ == '__main__':
    print(Run())
