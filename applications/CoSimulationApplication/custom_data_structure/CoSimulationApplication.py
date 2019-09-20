from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Custom "__init__.py" for the "CoSimulationApplication" module for the python-only version of the CoSimulationApplication

import os
kratos_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

co_simulation_py_scripts = os.path.join(kratos_path, "applications", "CoSimulationApplication", "python_scripts")

__path__.append(co_simulation_py_scripts) # adding the python-scripts from the CoSimulationApplication to the module