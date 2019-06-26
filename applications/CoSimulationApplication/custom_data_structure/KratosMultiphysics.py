from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

"""Custom "__init__.py" for the "KratosMultiphysics" module for the python-only version of the CoSimulationApplication
pyKratos is used to emulate the functionalities that are implemented in C++ in the Core
"""

import os
kratos_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

kratos_py_scripts = os.path.join(kratos_path, "kratos", "python_scripts")
pykratos = os.path.join(kratos_path, "applications", "CoSimulationApplication", "custom_data_structure", "pyKratos")

__path__.append(kratos_py_scripts) # adding the python-scripts from the core to the module
__path__.append(pykratos) # adding the python-scripts from pyKratos to the module

from .Parameters import Parameters
from .Model import Model
from .ModelPart import ModelPart
from .Variables import *
from .QuadElement import Quadrilateral3D4N
from .TriangleElement import Triangle
from .Logger import Logger

class KratosGlobals(object):
    def HasVariable(var_name):
        return False
    def GetVariable(var_name):
        return globals()[var_name]
    def __init__(self):
        pass

def Array1DVariable3(name):
    globals()[name+'_X'] = name+'_X'
    globals()[name+'_Y'] = name+'_Y'
    globals()[name+'_Z'] = name+'_Z'
    globals()[name] = [name, globals()[name+'_X'], globals()[name+'_Y'], globals()[name+'_Z']]
    return globals()[name]

def VariableDouble(name):
    globals()[name] = name