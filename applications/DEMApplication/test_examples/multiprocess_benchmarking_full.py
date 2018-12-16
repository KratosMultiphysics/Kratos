from __future__ import print_function
import os,subprocess,sys
import multiprocessing as mp
if sys.version_info >= (3, 0):
    import queue
else:
    import Queue as queue
from threading import Thread
import threading
from glob import glob
import shutil
import KratosMultiphysics.kratos_utilities as kratos_utils


kratos_benchmarking_path = os.path.join('..','..','..','benchmarking')
sys.path.append(kratos_benchmarking_path)
path = os.path.join('..','test_examples')
sys.path.append(path)
path = os.getcwd()
path = os.path.join(path,'basic_benchmarks')
os.chdir(path)
initial_number_of_threads = os.environ['OMP_NUM_THREADS']
os.environ['OMP_NUM_THREADS']='1'
os.system("echo Benchmarks will be running on $OMP_NUM_THREADS cpu")

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

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def run(benchmark):
    out_file_name = '{0}.info'.format(benchmark)
    f = open(out_file_name, 'wb')

    path_to_callable_script = os.path.join(path,"DEM_benchmarks.py")

    if sys.version_info >= (3, 0):
        subprocess.check_call(["python3", path_to_callable_script, str(benchmark)], stdout=f, stderr=f)
    else:
        subprocess.check_call(["python", "-3", path_to_callable_script, str(benchmark)], stdout=f, stderr=f)

    f.close()

    file_to_remove = out_file_name
    kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))

def worker(queue):
    """Process files from the queue."""
    for benchmark in iter(queue.get, None):
        try:
            print(Benchmark_text[benchmark - 1])
            run(benchmark)
        except Exception as e:# catch exceptions to avoid exiting the thread prematurely
            print("A problem was found in DEM Benchmark " + str(benchmark) + "... Resuming...\n")
            g = open("errors.err", "a")
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
            g.close()

def main():
    try:
        print("\nAdding processes to DEM parallel Benchmarking..............\n")
        g = open("errors.err", "w")
        g.write("The complete list of benchmarks are included at the end of this message as a quick reference.\n")
        g.write("\n========== DEM BENCHMARKING RESULTS ==========\n")
        g.write("\n=========== DEM DISCONTINUUM TESTS ===========\n")
        g.write("\n==== TSUJI PAPER BENCHMARKS. SLIDING REGIME ==\n\n")
        g.close()
        Text = ""
        failure = False

        q = queue.Queue()

        #Discontinuum Tests. From 1 to 17
        D_DEM_Benchmarks_list = list(range(1,18))

        #Continuum Tests
        C_DEM_Benchmarks_list = list(range(20,26))

        #Discontinuum Clusters Tests. From 30 to 33
        Dcl_DEM_Benchmarks_list = list(range(30,34))

        Total_DEM_Benchmarks_list = D_DEM_Benchmarks_list + C_DEM_Benchmarks_list + Dcl_DEM_Benchmarks_list

        for item in Total_DEM_Benchmarks_list:
            #print(Benchmark_text[item - 1])
            q.put_nowait(item)

        threads = [Thread(target=worker, args=(q,)) for _ in range(int(initial_number_of_threads))]
        for t in threads:
            t.daemon = True # threads die if the program dies
            t.start()
        for _ in threads: q.put_nowait(None) # signal no more files
        for t in threads: t.join() # wait for completion

        print('\n')
        g = open("errors.err", 'a')
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
        g.close()

        if 'FAILED' in open('errors.err').read():
            failure = True

        g = open("errors.err")
        file_contents = g.read()
        g.close()
        os.remove("errors.err")

        Text += file_contents.rstrip("\n")
        Text += "\n\n\n"
        delete_archives()
        return Text

    except:
        delete_archives()


def delete_archives():

    #.......................Removing extra files
    files_to_delete_list = glob('*.time')
    files_to_delete_list.extend(glob('*.dat'))
    files_to_delete_list.extend(glob('*.gp'))
    files_to_delete_list.extend(glob('*.txt'))
    files_to_delete_list.extend(glob('*.lst'))
    files_to_delete_list.extend(glob('*.info'))
    files_to_delete_list.extend(glob('*.err'))

    for to_erase_file in files_to_delete_list:
        os.remove(to_erase_file)

    #............Getting rid of unuseful folders
    folders_to_delete_list      = glob('*Data')
    folders_to_delete_list.extend(glob('*ists'))
    folders_to_delete_list.extend(glob('*ults'))
    folders_to_delete_list.extend(glob('*he__'))
    folders_to_delete_list.extend(glob('*aphs'))
    folders_to_delete_list.extend(glob('*iles'))

    for to_erase_folder in folders_to_delete_list:
        shutil.rmtree(to_erase_folder)




if __name__ == '__main__':
    mp.freeze_support() # optional if the program is not frozen
    print(main())
