from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys
import platform

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
swimming_dem_scripts_path = 'hydrodynamic_forces'
sys.path.append(swimming_dem_scripts_path)
import benchmarking
os.chdir(swimming_dem_scripts_path)


def Run():
    
    print("\nStarting swimming_DEM Benchmarks..............\n")    
    Text=""           
        
    if platform.system()=="Windows":
        os.system("python hydrodynamic_forces.py " + " > BenchTemp.txt")
    else:
        if sys.version_info >= (3, 0):
            os.system("python3 hydrodynamic_forces.py " + " > BenchTemp.txt")
        else:
            os.system("python -3 hydrodynamic_forces.py " + " > BenchTemp.txt")
                
    os.remove("BenchTemp.txt")
        
    f = open("hydrodynamic_forces.txt")
    file_contents = f.read()
    f.close()
    
    Text += file_contents.rstrip("\n")
    Text += "\n\n\n"        
    
    return Text
       
if __name__ == '__main__':
    print(Run())
