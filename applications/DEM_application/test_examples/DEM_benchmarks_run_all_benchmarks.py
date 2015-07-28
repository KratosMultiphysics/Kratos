from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys
import platform

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
dem_scripts_path = '../test_examples/basic_benchmarks'
sys.path.append(dem_scripts_path)
import basic_benchmarks
import benchmarking
os.chdir(dem_scripts_path)


def Run():
    
    Text  = "\n========== DEM BENCHMARKS ===========\n"
    Text += "========== SLIDING REGIME ===========\n\n"
    
    for benchmark in range(3, 9):
        
        error1, error2, error3 = basic_benchmarks.Run(benchmark)
        
        Text += "Test " + str(benchmark) + ":"

        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            Text += " OK!........ Test " + str(benchmark) + " SUCCESSFUL\n"
        else:
            Text += " KO!........ Test " + str(benchmark) + " FAILED\n"
    
    Text +="\n\n"
    
    return Text

    
def RunBenchmark(ExamplePath):

    if platform.system()=="Windows":
        os.system("python " + ExamplePath)
    else:
        if sys.version_info >= (3, 0):
            os.system(
                "python3 " +
                ExamplePath)
        else:
            os.system(
                "python -3 " +
                ExamplePath)
       
if __name__ == '__main__':
    print(Run())
