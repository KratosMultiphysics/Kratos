import os
import sys
import platform
from KratosMultiphysics.testing.utilities import GetPython3Command

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
swimming_dem_scripts_path = 'hydrodynamic_forces'
sys.path.append(swimming_dem_scripts_path)
import benchmarking
os.chdir(swimming_dem_scripts_path)


def Run():

    print("\nStarting swimming_DEM Benchmarks..............\n")
    Text=""

    py_cmd = GetPython3Command()
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
