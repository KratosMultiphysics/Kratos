import sys
import os

# sys.path.append(os.path.join('..', '..', '..'))
# sys.path.append(os.path.join('..', 'python_scripts'))
# sys.path.append(r"D:\software_development\Kratos\bin\Debug")

from KratosMultiphysics import Tester
import KratosMultiphysics.GeoMechanicsApplication

# If you want to run tests defined in an application, import it here.

Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS) # Set the verbosity level

#Tester.RunAllTestCases() #Test all cases

Tester.RunTestSuite("KratosGeoMechanicsFastSuite") #Test a whole suite

#Tester.RunTestCases("Test_name_here") #Test a specific case