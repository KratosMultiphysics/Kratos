# pyKratos imports
from .Variable import CreateDoubleVariable

def RegisterVariables(module):
    CreateDoubleVariable(module, "SCALAR_DISPLACEMENT")
    CreateDoubleVariable(module, "SCALAR_ROOT_POINT_DISPLACEMENT")
    CreateDoubleVariable(module, "SCALAR_REACTION")
    CreateDoubleVariable(module, "SCALAR_FORCE")
    CreateDoubleVariable(module, "SCALAR_VOLUME_ACCELERATION")
