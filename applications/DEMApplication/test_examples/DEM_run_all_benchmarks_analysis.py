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

Benchmark_text = ["Running DEM Benchmark 1.... Elastic normal impact of two identical spheres\n",
                  "Running DEM Benchmark 2.... Elastic normal impact of a sphere against a rigid plane\n",
                  "Running DEM Benchmark 3.... Impact of a sphere against a rigid plane with different coefficients of restitution\n",
                  "Running DEM Benchmark 4.... Oblique impact of a sphere with a rigid plane with constant velocity module and variable incident angles\n",
                  "Running DEM Benchmark 5.... Oblique impact of a sphere with a rigid plane with constant normal velocity and different tangential velocities\n",
                  "Running DEM Benchmark 6.... Impact of a sphere with a rigid plane with a constant normal velocity and variable angular velocities\n",
                  "Running DEM Benchmark 7.... Impact of two identical spheres with a constant normal velocity and different angular velocities\n",
                  "Running DEM Benchmark 8.... Impact of two differently sized spheres with a constant normal velocity and variable angular velocities\n",
                  "Running DEM Benchmark 9.... Impact of two identical spheres with a constant normal velocity and different coefficients of restitution\n",
                  "Running DEM Benchmark 10... Linear: Oblique impact of a sphere with an elastic plane with constant normal velocity and different angular velocities\n",
                  "Running DEM Benchmark 11... Hertzian: Oblique impact of a sphere with an elastic plane with constant normal velocity and different angular velocities\n",
                  "Running DEM Benchmark 12... Sphere rotating over a plane surface with Rolling Friction\n",
                  "Running DEM Benchmark 13... Impact of a low stiffness sphere against a rigid plane divided in small triangular elements\n",
                  "Running DEM Benchmark 14... Impact of a low stiffness sphere against a rigid edge divided in small triangular elements\n",
                  "Running DEM Benchmark 15... Impact of a low stiffness sphere against a rigid vertex divided in small triangular elements\n",
                  "Running DEM Benchmark 16... Spheres contacting multiple entities (facets, edges and vertices)\n",
                  "Running DEM Benchmark 17... Sphere sliding on a plane (discretized with triangles and quadrilaterals) with friction\n",
                  "","",
                  "Running DEM Benchmark 20... Normal compression of two identical spheres\n",\
                  "Running DEM Benchmark 21... Normal compression of two identical indented spheres\n",\
                  "Running DEM Benchmark 22... Tensile test of two identical spheres\n",\
                  "Running DEM Benchmark 23... Tensile test of two identical indented spheres\n",\
                  "Running DEM Benchmark 24... Shear test of two identical spheres by applying rotation\n",\
                  "Running DEM Benchmark 25... Shear test of two identical spheres by applying rotation and radius expansion\n",\
                  "","","","",
                  "Running DEM Benchmark 30... Cylinder cluster with imposed angular velocity in two axis (Velocity Verlet + Zhao scheme)\n",
                  "Running DEM Benchmark 31... Cylinder cluster with imposed angular velocity in two axis (Symplectic Euler + Runge-Kutta scheme)\n",
                  "Running DEM Benchmark 32... Fiber cluster bouncing without any damping (Velocity Verlet + Zhao scheme)\n",
                  "Running DEM Benchmark 33... Fiber cluster bouncing without any damping (Symplectic Euler + Runge-Kutta scheme)\n"]

def Run():

    print("\nStarting DEM Benchmarking..............\n")

    with open("errors.err", "w") as g:
        g.write("The complete list of benchmarks are included at the end of this message as a quick reference.\n")
        g.write("\n========== DEM BENCHMARKING RESULTS ==========\n")
        g.write("\n=========== DEM DISCONTINUUM TESTS ===========\n")
        g.write("\n==== TSUJI PAPER BENCHMARKS. SLIDING REGIME ==\n\n")
    Text = ""

    #list_of_failed_tests = []

    #Discontinuum Tests. From 1 to 17
    D_DEM_Benchmarks_list = list(range(1,18))
    #Continuum Tests
    C_DEM_Benchmarks_list = list(range(20,26))
    #Discontinuum Clusters Tests. From 30 to 33
    Dcl_DEM_Benchmarks_list = list(range(30,34))

    Total_DEM_Benchmarks_list = D_DEM_Benchmarks_list + C_DEM_Benchmarks_list + Dcl_DEM_Benchmarks_list

    failure = False

    with open("BenchTemp.info", "w") as f:
        for benchmark in Total_DEM_Benchmarks_list:
            print(Benchmark_text[benchmark - 1])
            try:
                os.environ["OMP_NUM_THREADS"] = "1"

                if platform.system()=="Windows":
                    subprocess.check_call(["python", path + "/DEM_benchmarks_analysis.py", str(benchmark), ">", "BenchTemp.info"], stdout=f, stderr=f)
                else:
                    if sys.version_info >= (3, 0):
                        subprocess.check_call(["python3", path + "/DEM_benchmarks_analysis.py", str(benchmark), ">", "BenchTemp.info"], stdout=f, stderr=f)
                    else:
                        subprocess.check_call(["python", "-3", path + "/DEM_benchmarks_analysis.py", str(benchmark), ">", "BenchTemp.info"], stdout=f, stderr=f)

                os.environ["OMP_NUM_THREADS"] = "" # Trying to set a 'default' value
            except:
                #failure = True
                #list_of_failed_tests += [benchmark]
                print("A problem was found in DEM Benchmark " + str(benchmark) + "... Resuming...\n")
                with open("errors.err", "a") as g:
                    if benchmark == 10:
                        g.write("\n===== THORNTON PAPER TESTS. FULL REGIME. LINEAR LAW =====\n\n")
                    if benchmark == 11:
                        g.write("\n===== THORNTON PAPER TESTS. FULL REGIME. HERTZIAN LAW ===\n\n")
                    if benchmark == 12:
                        g.write("\n===== WENSRICH PAPER TEST. ROLLING FRICTION =============\n\n")
                    if benchmark == 13:
                        g.write("\n===== DE/FE CONTACT BENCHMARKS ==========================\n\n")
                    if benchmark == 20:
                        g.write("\n===== BASIC CONTINUUM TESTS  ============================\n\n")
                    if benchmark == 30:
                        g.write("\n===== DISCONTINUUM CLUSTERS TESTS  ======================\n\n")
                    g.write("DEM Benchmark " + str(benchmark) + ": KO!........ Test " + str(benchmark) + " FAILED\n")
        print('\n')


    with open("errors.err", "a") as g:
        g.write("\n---------------------------------------------------------------------\n")
        g.write("\nList of Benchmarks:\n")
        g.write("\nDISCONTINUUM TESTS:\n")
        g.write("Benchmark 01. Elastic normal impact of two identical spheres\n")
        g.write("Benchmark 02. Elastic normal impact of a sphere against a rigid plane\n")
        g.write("Benchmark 03. Impact of a sphere against a rigid plane with different coefficients of restitution\n")
        g.write("Benchmark 04. Oblique impact of a sphere with a rigid plane with constant velocity module and variable incident angles\n")
        g.write("Benchmark 05. Oblique impact of a sphere with a rigid plane with constant normal velocity and different tangential velocities\n")
        g.write("Benchmark 06. Impact of a sphere with a rigid plane with a constant normal velocity and variable angular velocities\n")
        g.write("Benchmark 07. Impact of two identical spheres with a constant normal velocity and different angular velocities\n")
        g.write("Benchmark 08. Impact of two differently sized spheres with a constant normal velocity and variable angular velocities\n")
        g.write("Benchmark 09. Impact of two identical spheres with a constant normal velocity and different coefficients of restitution\n")
        g.write("Benchmark 10. Oblique impact of a sphere with an elastic plane with constant normal velocity and different angular velocities\n")
        g.write("Benchmark 11. Oblique impact of a sphere with an elastic plane with constant normal velocity and different angular velocities\n")
        g.write("Benchmark 12. Sphere rotating over a plane surface with Rolling Friction\n")
        g.write("Benchmark 13. Impact of a low stiffness sphere against a rigid plane divided in small triangular elements\n")
        g.write("Benchmark 14. Impact of a low stiffness sphere against a rigid edge divided in small triangular elements\n")
        g.write("Benchmark 15. Impact of a low stiffness sphere against a rigid vertex divided in small triangular elements\n")
        g.write("Benchmark 16. Spheres contacting multiple entities (facets, edges and vertices)\n")
        g.write("Benchmark 17. Sphere sliding on a plane (discretized with triangles and quadrilaterals) with friction\n")
        g.write("\nCONTINUUM TESTS:\n")
        g.write("Benchmark 20. Normal compression of two identical spheres\n")
        g.write("Benchmark 21. Normal compression of two identical indented spheres\n")
        g.write("Benchmark 22. Tensile test of two identical spheres\n")
        g.write("Benchmark 23. Tensile test of two identical indented spheres\n")
        g.write("Benchmark 24. Shear test of two identical spheres by applying rotation\n")
        g.write("Benchmark 25. Shear test of two identical spheres by applying rotation and radius expansion\n")
        g.write("\nDISCONTINUUM CLUSTERS TESTS:\n")
        g.write("Benchmark 30. Cylinder cluster with imposed angular velocity in two axis (Velocity Verlet + Zhao scheme)\n")
        g.write("Benchmark 31. Cylinder cluster with imposed angular velocity in two axis (Symplectic Euler + Runge-Kutta scheme)\n")
        g.write("Benchmark 32. Fiber cluster bouncing without any damping (Velocity Verlet + Zhao scheme)\n")
        g.write("Benchmark 33. Fiber cluster bouncing without any damping (Symplectic Euler + Runge-Kutta scheme)\n")

    with open("errors.err") as g:
        file_contents = g.read()
        if 'FAILED' in file_contents:
            failure = True

    os.remove("errors.err")

    Text += file_contents.rstrip("\n")
    Text += "\n\n\n"

    # To send an email summary to the DEM Team
    #list_of_failed_tests = ', '.join(str(e) for e in list_of_failed_tests)
    #subject = "DEM Benchmarks Results. Problems found in tests number "
    #subject += list_of_failed_tests
    recipients = ["latorre@cimne.upc.edu"]
    recipients += ["maceli@cimne.upc.edu", "gcasas@cimne.upc.edu", "farrufat@cimne.upc.edu", "jirazabal@cimne.upc.edu"]
    subject = "Problems found in DEM Benchmarks"
    message = "From: Kratos Benchmarking <no-reply-kratos-benchmarking@cimne.upc.es>\nSubject: " + subject + "\n" + Text

    #if failure:
    #if (__name__ != '__main__' and failure):
    if (False):
        smtplib.SMTP("smtps.cimne.upc.es").sendmail("Kratos Benchmarking <no-reply-kratos-benchmarking@cimne.upc.es>", recipients, message)
    return Text


if __name__ == '__main__':
    print(Run())
