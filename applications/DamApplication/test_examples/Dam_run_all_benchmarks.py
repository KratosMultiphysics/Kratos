from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import subprocess
import sys
import platform
import smtplib

path = '../test_examples'
sys.path.append(path)

path = os.getcwd()
path += '/basic_benchmarks'
os.chdir(path)

Benchmark_text = ["Running Dam Benchmark 2D 01... Mechanical-PlaneStrain with hydrostatic pressure\n",
                  "Running Dam Benchmark 2D 02... Mechanical-PlaneStress with hydrostatic pressure\n",
                  "Running Dam Benchmark 2D 03... ThermoMechanical-PlaneStrain with hydrostatic pressure and heat flux\n",
                  "Running Dam Benchmark 2D 04... ThermoMechanical-PlaneStress with hydrostatic pressure and heat flux\n",
                  "Running Dam Benchmark 2D 05... ThermoMechanical-PlaneStrain with hydrostatic pressure and imposed temperature\n",
                  "Running Dam Benchmark 2D 06... ThermoMechanical-PlaneStress with hydrostatic pressure and imposed temperature\n",
                  "Running Dam Benchmark 2D 07... ThermoMechanical-PlaneStrain with hydrostatic pressure and Bofang\n",
                  "Running Dam Benchmark 2D 08... ThermoMechanical-PlaneStress with hydrostatic pressure and Bofang\n",
                  "Running Dam Benchmark 2D 09... ThermoMechanical-PlaneStrain with heat flux\n",
                  "Running Dam Benchmark 2D 10... ThermoMechanical-PlaneStress with heat flux\n",
                  "Running Dam Benchmark 2D 11... ThermoMechanical-PlaneStrain with imposed temperature\n",
                  "Running Dam Benchmark 2D 12... ThermoMechanical-PlaneStress with imposed temperature\n\n",
                  "Running Dam Benchmark 3D 01... Mechanical-PlaneStrain with hydrostatic pressure\n",\
                  "Running Dam Benchmark 3D 02... ThermoMechanical with hydrostatic pressure and heat flux\n",\
                  "Running Dam Benchmark 3D 03... ThermoMechanical with hydrostatic pressure and imposed temperature\n",\
                  "Running Dam Benchmark 3D 04... ThermoMechanical with hydrostatic pressure and Bofang\n",\
                  "Running Dam Benchmark 3D 05... ThermoMechanical with heat flux\n",\
                  "Running Dam Benchmark 3D 06... ThermoMechanical with imposed temperature\n"]

def Run():

    print("\nStarting Dam Benchmarking..............\n")

    with open("errors.err", "w") as g:
        g.write("The complete list of benchmarks are included at the end of this message as a quick reference.\n")
        g.write("\n========== DAM BENCHMARKING RESULTS ==========\n\n")
    Text = ""

    Benchmarks_list_2D = list(range(201,213))
    Benchmarks_list_3D = list(range(301,307))

    Total_Dam_Benchmarks_list = Benchmarks_list_2D + Benchmarks_list_3D

    failure = False

    with open("BenchTemp.info", "w") as f:
        i = 0
        for benchmark in Total_Dam_Benchmarks_list:
            print(Benchmark_text[i])
            i += 1
            try:
                if platform.system()=="Windows":
                    os.system("setenv OMP_NUM_THREADS 1") # Is that the correct way to run on Windows?
                    subprocess.check_call(["python", path + "/Dam_benchmarks.py", str(benchmark), ">", "BenchTemp.info"], stdout=f, stderr=f)
                    os.system("setenv OMP_NUM_THREADS ") # Trying to set a 'default' value

                else:
                    os.environ['OMP_NUM_THREADS']='1'
                    if sys.version_info >= (3, 0):
                        subprocess.check_call(["python3", path + "/Dam_benchmarks.py", str(benchmark), ">", "BenchTemp.info"], stdout=f, stderr=f)

                    else:
                        subprocess.check_call(["python", "-3", path + "/Dam_benchmarks.py", str(benchmark), ">", "BenchTemp.info"], stdout=f, stderr=f)
                    os.system("export OMP_NUM_THREADS=") # Trying to set a 'default' value
            except:
                #failure = True
                #list_of_failed_tests += [benchmark]
                print("A problem was found in DEM Benchmark " + str(benchmark) + "... Resuming...\n")
                with open("errors.err", "a") as g:
                    if benchmark == 201:
                        g.write("\n===================== 2D BENCHMARKS =====================\n\n")
                    if benchmark == 301:
                        g.write("\n===================== 3D BENCHMARKS =====================\n\n")
                    g.write("DEM Benchmark " + str(benchmark) + ": KO!........ Test " + str(benchmark) + " FAILED\n")
        print('\n')

    with open("errors.err", "a") as g:
        g.write("\n---------------------------------------------------------------------\n")
        g.write("\nList of Benchmarks:\n")
        g.write("\nTWO-DIMENSIONAL TESTS:\n")
        g.write("Dam Benchmark 2D 01: Mechanical-PlaneStrain with hydrostatic pressure\n")
        g.write("Dam Benchmark 2D 02: Mechanical-PlaneStress with hydrostatic pressure\n")
        g.write("Dam Benchmark 2D 03: ThermoMechanical-PlaneStrain with hydrostatic pressure and heat flux\n")
        g.write("Dam Benchmark 2D 04: ThermoMechanical-PlaneStress with hydrostatic pressure and heat flux\n")
        g.write("Dam Benchmark 2D 05: ThermoMechanical-PlaneStrain with hydrostatic pressure and imposed temperature\n")
        g.write("Dam Benchmark 2D 06: ThermoMechanical-PlaneStress with hydrostatic pressure and imposed temperature\n")
        g.write("Dam Benchmark 2D 07: ThermoMechanical-PlaneStrain with hydrostatic pressure and Bofang\n")
        g.write("Dam Benchmark 2D 08: ThermoMechanical-PlaneStress with hydrostatic pressure and Bofang\n")
        g.write("Dam Benchmark 2D 09: ThermoMechanical-PlaneStrain with heat flux\n")
        g.write("Dam Benchmark 2D 10: ThermoMechanical-PlaneStress with heat flux\n")
        g.write("Dam Benchmark 2D 11: ThermoMechanical-PlaneStrain with imposed temperature\n")
        g.write("Dam Benchmark 2D 12: ThermoMechanical-PlaneStress with imposed temperature\n")
        g.write("\nTHREE-DIMENSIONAL TESTS:\n")
        g.write("Dam Benchmark 3D 01: Mechanical-PlaneStrain with hydrostatic pressure\n")
        g.write("Dam Benchmark 3D 02: ThermoMechanical with hydrostatic pressure and heat flux\n")
        g.write("Dam Benchmark 3D 03: ThermoMechanical with hydrostatic pressure and imposed temperature\n")
        g.write("Dam Benchmark 3D 04: ThermoMechanical with hydrostatic pressure and Bofang\n")
        g.write("Dam Benchmark 3D 05: ThermoMechanical with heat flux\n")
        g.write("Dam Benchmark 3D 06: ThermoMechanical with imposed temperature\n")

    with open("errors.err") as g:
        file_contents = g.read()
        if 'FAILED' in file_contents:
            failure = True

    os.remove("errors.err")

    Text += file_contents.rstrip("\n")
    Text += "\n\n\n"
    
    return Text

if __name__ == '__main__':
    print(Run())
