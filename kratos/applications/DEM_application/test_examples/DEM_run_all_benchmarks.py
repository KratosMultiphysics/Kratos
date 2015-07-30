from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import subprocess
import sys
import platform

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
path = '../test_examples'
sys.path.append(path)
import benchmarking
path = os.getcwd()
path += '/basic_benchmarks'
os.chdir(path)


def Run():
    
    print("\nStarting DEM Benchmarks..............\n")
    g = open("errors.txt", "w")
    g.write("\n========== DEM BENCHMARKS ===========\n")
    g.write("========== SLIDING REGIME ===========\n\n")
    g.close()
    Text=""
    
    f=open("BenchTemp.txt", "w")

    for benchmark in range(3, 6):
      
        print("Running Benchmark " + str(benchmark) + " of 8.............")
             
        try:
            if platform.system()=="Windows":
                subprocess.check_call(["python", path + "/Chung_Ooi_benchmarks.py", str(benchmark), ">", "BenchTemp.txt"], stdout=f, stderr=f)
                            
            else:
                if sys.version_info >= (3, 0):
                    subprocess.check_call(["python3", path + "/Chung_Ooi_benchmarks.py", str(benchmark), ">", "BenchTemp.txt"], stdout=f, stderr=f)
                    
                else:
                    subprocess.check_call(["python", "-3", path + "/Chung_Ooi_benchmarks.py", str(benchmark), ">", "BenchTemp.txt"], stdout=f, stderr=f)
        except:
            print("\nDEM Benchmarks: A problem was found in Benchmark " + str(benchmark) + "... Resuming...\n")
    
    print('\n')
    f.close()
    os.remove("BenchTemp.txt")        
    
    g = open("errors.txt")
    file_contents = g.read()
    g.close()
    os.remove("errors.txt")
    
    Text += file_contents.rstrip("\n")
    Text += "\n\n\n"
        
    return Text


if __name__ == '__main__':
    print(Run())
