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
                  
Benchmark_text = ["Running DEM Benchmark 1... Elastic normal impact of two identical spheres\n",
                  "Running DEM Benchmark 2... Elastic normal impact of a sphere against a rigid plane\n",
                  "Running DEM Benchmark 3... Impact of a sphere against a rigid plane with different coefficients of restitution\n",
                  "Running DEM Benchmark 4... Oblique impact of a sphere with a rigid plane with constant velocity module and variable incident angles\n",
                  "Running DEM Benchmark 5... Oblique impact of a sphere with a rigid plane with constant normal velocity and different angular velocities\n",
                  "Running DEM Benchmark 6... Impact of a sphere with a rigid plane with a constant normal velocity and variable angular velocities\n",
                  "Running DEM Benchmark 7... Impact of two identical spheres with a constant normal velocity and different angular velocities\n",
                  "Running DEM Benchmark 8... Impact of two differently sized spheres with a constant normal velocity and variable angular velocities\n",
                  "Running DEM Benchmark 9... Impact of two identical spheres with a constant normal velocity and different coefficients of restitution\n",
                  "","","","","","","","","","",
                  "Running DEM Benchmark 20... Normal compression of two identical spheres\n",\
                  "Running DEM Benchmark 21... Normal compression of two identical indented spheres\n",\
                  "Running DEM Benchmark 22... Tensile\n",\
                  "Running DEM Benchmark 23... Indented tensile\n",\
                  "Running DEM Benchmark 24... Shear\n",\
                  "Running DEM Benchmark 25... Shear + radius expansion\n"]

def Run():
    
    print("\nStarting DEM Benchmarking..............\n")
    g = open("errors.txt", "w")
    g.write("\n========== DEM BENCHMARKING RESULTS ==========\n\n")
    g.write("== BASIC DISCONTINUUM TESTS, SLIDING REGIME ==\n\n")
    g.close()
    Text=""
    f=open("BenchTemp.txt", "w")
    
    #Discontinuum Tests. From 1 to 9
    D_DEM_Benchmarks_list = list(range(1,10))

    #Continuum Tests
    C_DEM_Benchmarks_list = list(range(20,26))
    
    Total_DEM_Benchmarks_list =  D_DEM_Benchmarks_list + C_DEM_Benchmarks_list
    
    for benchmark in Total_DEM_Benchmarks_list:
          
        print(Benchmark_text[benchmark - 1])
                
        try:
            if platform.system()=="Windows":
                os.system("setenv OMP_NUM_THREADS 1") # Is that the correct way to run on Windows?
                subprocess.check_call(["python", path + "/DEM_benchmarks.py", str(benchmark), ">", "BenchTemp.txt"], stdout=f, stderr=f)
                os.system("setenv OMP_NUM_THREADS 16") # Trying to set a 'default' value
                
            else:
                os.system("export OMP_NUM_THREADS=1")
                if sys.version_info >= (3, 0):
                    subprocess.check_call(["python3", path + "/DEM_benchmarks.py", str(benchmark), ">", "BenchTemp.txt"], stdout=f, stderr=f)
                    
                else:
                    subprocess.check_call(["python", "-3", path + "/DEM_benchmarks.py", str(benchmark), ">", "BenchTemp.txt"], stdout=f, stderr=f)
                os.system("export OMP_NUM_THREADS=16") # Trying to set a 'default' value
        except:
            print("A problem was found in DEM Benchmark " + str(benchmark) + "... Resuming...\n")
            g = open("errors.txt", "a")
            g.write("DEM Benchmark " + str(benchmark) + ": KO!........ Test " + str(benchmark) + " FAILED\n")
            g.close()
            
    print('\n')
    f.close()
    os.remove("BenchTemp.txt")        

    g = open("errors.txt")
    file_contents = g.read()
    g.close()
    os.remove("errors.txt")
    
    Text += file_contents.rstrip("\n")
    Text += "\n\n\n"

    '''
    print("\nStarting DEM Benchmarking..............\n")
    g2 = open("errors2.txt", "w")
    g2.write("\n========== DEM BENCHMARKING RESULTS ==========\n\n")
    g2.write("== BASIC CONTINUUM TESTS ==\n\n")
    g2.close()
    Text2=""
    f2=open("BenchTemp2.txt", "w")




    Benchmark_text2 = ["Running DEM Benchmark 20... Normal compression of two identical spheres\n",\
                      "Running DEM Benchmark 21... Normal compression of two identical indented spheres\n",\
                      "Running DEM Benchmark 22... Tensile\n",\
                      "Running DEM Benchmark 23... Indented tensile\n",\
                      "Running DEM Benchmark 24... Shear\n",\
                      "Running DEM Benchmark 25... Shear + radius expansion\n"]

    for benchmark in range(20, 25):

        print(Benchmark_text2[benchmark - 20])

        try:
            if platform.system()=="Windows":
                subprocess.check_call(["python", path + "/DEM_benchmarks.py", str(benchmark), ">", "BenchTemp2.txt"], stdout=f, stderr=f)

            else:
                if sys.version_info >= (3, 0):
                    subprocess.check_call(["python3", path + "/DEM_benchmarks.py", str(benchmark), ">", "BenchTemp2.txt"], stdout=f, stderr=f)

                else:
                    subprocess.check_call(["python", "-3", path + "/DEM_benchmarks.py", str(benchmark), ">", "BenchTemp2.txt"], stdout=f, stderr=f)
        except:
            print("A problem was found in DEM Benchmark " + str(benchmark) + "... Resuming...\n")
            g2 = open("errors2.txt", "a")
            g2.write("DEM Benchmark " + str(benchmark) + ": KO!........ Test " + str(benchmark) + " FAILED\n")
            g2.close()


    print('\n')
    f2.close()
    os.remove("BenchTemp2.txt")

    g2 = open("errors2.txt")
    file_contents2 = g2.read()
    g2.close()
    os.remove("errors2.txt")

    Text2 += file_contents2.rstrip("\n")
    Text2 += "\n\n\n"
    '''
    return Text


if __name__ == '__main__':
    print(Run())
