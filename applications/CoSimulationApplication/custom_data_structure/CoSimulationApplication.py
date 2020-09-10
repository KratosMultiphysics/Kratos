from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Custom "__init__.py" for the "CoSimulationApplication" module for the python-only version of the CoSimulationApplication

import os, sys
kratos_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

co_simulation_py_scripts = os.path.join(kratos_path, "applications", "CoSimulationApplication", "python_scripts")

__path__.append(co_simulation_py_scripts) # adding the python-scripts from the CoSimulationApplication to the module

from ..variables_co_simulation import RegisterVariables

print("""
    KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ __
           | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_ \\
           | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
            \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
""")

RegisterVariables(sys.modules[__name__])
