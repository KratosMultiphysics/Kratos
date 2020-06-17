from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys
import platform
import shutil

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
swimming_dem_scripts_path = 'hydrodynamic_forces'
sys.path.append(swimming_dem_scripts_path)
import benchmarking
os.chdir(swimming_dem_scripts_path)


def Run():

    print("\nStarting swimming_DEM Benchmarks..............\n")
    Text=""

    py_cmd = "python3" if shutil.which("python3") is not None else "python"
    os.system(py_cmd + " hydrodynamic_forces.py " + " > BenchTemp.txt")

    os.remove("BenchTemp.txt")

    f = open("hydrodynamic_forces.txt")
    file_contents = f.read()
    f.close()

    Text += file_contents.rstrip("\n")
    Text += "\n\n\n"

    return Text

if __name__ == '__main__':
    print(Run())
