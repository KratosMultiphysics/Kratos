from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Custom "__init__.py" for the "KratosMultiphysics" module for the python-only version of the CoSimulationApplication
# pyKratos is used to emulate the functionalities that are implemented in C++ in the Core

import os, sys
kratos_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

kratos_py_scripts = os.path.join(kratos_path, "kratos", "python_scripts")
pykratos = os.path.join(kratos_path, "applications", "CoSimulationApplication", "custom_data_structure", "pyKratos")

__path__.append(kratos_py_scripts) # adding the python-scripts from the core to the module
__path__.append(pykratos) # adding the python-scripts from pyKratos to the module

# bringing all classes to the scope of the Module
# this way e.g. "from KratosMultiphysics import ModelPart" becomes possible
# => this way pyKratos is fully compatible with the real Kratos
from .Parameters import Parameters
from .Model import Model
from .ModelPart import ModelPart
from .variables_kratos import RegisterVariables
from .QuadriateralElement import QuadrilateralElement
from .TriangleElement import TriangleElement
from .Logger import Logger
from .IntervalUtility import IntervalUtility

print("""              _  __          _
  _ __  _   _| |/ /_ __ __ _| |_ ___  ___
 | '_ \| | | | ' /| '__/ _` | __/ _ \/ __|
 | |_) | |_| | . \| | | (_| | || (_) \__ \\
 | .__/ \__, |_|\_\_|  \__,_|\__\___/|___/
 |_|    |___/
""", end='')

_registered_variables = {}

class KratosGlobals(object):
    def HasVariable(var_name):
        return var_name in _registered_variables

    def GetVariable(var_name):
        if not KratosGlobals.HasVariable(var_name):
            raise Exception('Variable "{}" does not exist!'.format(var_name))
        return _registered_variables[var_name]

    def GetVariableType(var_name):
        if not KratosGlobals.HasVariable(var_name):
            return "None"
        return KratosGlobals.GetVariable(var_name).Type()

def IsDistributedRun():
    return False # pyKratos cannot be run in MPI

RegisterVariables(sys.modules[__name__])
