from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import subprocess
import sys
import platform
import smtplib

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
path = '../test_examples'
sys.path.append(path)
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
                  "Running DEM Benchmark 10... Linear: Oblique impact of a sphere with an elastic plane with constant normal velocity and different angular velocities\n",
                  "Running DEM Benchmark 11... Hertzian: Oblique impact of a sphere with an elastic plane with constant normal velocity and different angular velocities\n",
                  "","","","","","","","",
                  "Running DEM Benchmark 20... Normal compression of two identical spheres\n",\
                  "Running DEM Benchmark 21... Normal compression of two identical indented spheres\n",\
                  "Running DEM Benchmark 22... Tensile test of two identical spheres\n",\
                  "Running DEM Benchmark 23... Tensile test of two identical indented spheres\n",\
                  "Running DEM Benchmark 24... Shear test of two identical spheres by applying rotation\n",\
                  "Running DEM Benchmark 25... Shear test of two identical spheres by applying rotation and radius expansion\n"]

def Run():
    
    print("\nStarting DEM Benchmarking..............\n")
    g = open("errors.txt", "w")
    g.write("\n\n========== DEM BENCHMARKING RESULTS ==========\n\n")
    g.write("\n=========== DEM DISCONTINUUM TESTS ===========\n")
    g.write("\n==== TSUJI PAPER BENCHMARKS. SLIDING REGIME ==\n\n")
    g.close()
    Text = ""
    f = open("BenchTemp.txt", "w")
    failure = False
    
    #Discontinuum Tests. From 1 to 11
    D_DEM_Benchmarks_list = list(range(1,12))
        
    #Continuum Tests
    C_DEM_Benchmarks_list = list(range(20,26))
        
    Total_DEM_Benchmarks_list = D_DEM_Benchmarks_list + C_DEM_Benchmarks_list
    
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
            failure = True
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
    
    # To send an email summary to the DEM Team
    recipients = ["latorre@cimne.upc.edu"]
    recipients += ["maceli@cimne.upc.edu", "msantasusana@cimne.upc.edu", "gcasas@cimne.upc.edu", "farrufat@cimne.upc.edu", "jirazabal@cimne.upc.edu"]
    subject = "Problems found in DEM Benchmarks"   
    message = "From: Kratos Benchmarking <no-reply-kratos-benchmarking@cimne.upc.es>\nSubject: " + subject + "\n" + Text
    
    if (__name__ != '__main__' and failure):
        smtplib.SMTP("smtps.cimne.upc.es").sendmail("Kratos Benchmarking <no-reply-kratos-benchmarking@cimne.upc.es>", recipients, message)
    
    return Text


if __name__ == '__main__':
    print(Run())
