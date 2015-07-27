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

def Run():
    
    Text = "\n========== DEM BENCHMARKS ==========\n\n"
    error1_3, error2_3, error3_3 = basic_benchmarks.Run(3)
    Text += "Test 3: "
    
    if (error1_3 < 10.0):
        Text += "OK!........"
        Text += "Test 3 SUCCESSFUL\n"
    else:
        Text += "KO!............"
        Text += "Test 3 FAILED\n"

    error1_4, error2_4, error3_4 = basic_benchmarks.Run(4)
    Text += "Test 4: "
    
    if (error1_4 < 10.0 and error2_4 < 10.0 and error3_4 < 10.0):
        Text += "OK!........"
        Text += "Test 4 SUCCESSFUL\n"
    else:
        Text += "KO!............"
        Text += "Test 4 FAILED\n"

    error1_7, error2_7, error3_7 = basic_benchmarks.Run(7)
    Text += "Test 7: "
    
    if (error1_7 < 10.0 and error2_7 < 10.0):
        Text += "OK!........"
        Text += "Test 7 SUCCESSFUL\n\n\n"
    else:
        Text += "KO!............"
        Text += "Test 7 FAILED\n\n\n"
    
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
    Run()
