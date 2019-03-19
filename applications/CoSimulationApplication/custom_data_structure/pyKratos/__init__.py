from .Parameters import Parameters
from .Model import Model
from .ModelPart import ModelPart
#from .Vector import Vector
#from .Matrix import Matrix
from .Variables import *
from .QuadElement import Quadrilateral3D4N
from .TriangleElement import Triangle

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