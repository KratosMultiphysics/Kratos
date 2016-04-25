from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "========== Kratos Core ==========\n"

    #
    # Check if all the variables are registered in python -- only issues a warning
    Text += "Registered In Python test: "
    #os.chdir(".")
    sys.path.append(os.getcwd())
    
    import test_variables_python_interface
    Text += test_variables_python_interface.CheckVariables()
    #os.chdir("..")

    
    ##check geometries
    #import KratosMultiphysics
    #Text += KratosMultiphysics.GeometryTesterUtility().RunTest()
    
    
    
    ###### here final output
    print(Text)
    return Text

if __name__ == '__main__':
    Run()
