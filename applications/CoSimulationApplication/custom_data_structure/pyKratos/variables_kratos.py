from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# pyKratos imports
from .Variable import CreateArray3Variable
from .Variable import CreateDoubleVariable
from .Variable import CreateVectorVariable

def RegisterVariables(module):
    CreateArray3Variable(module, "DISPLACEMENT")
    CreateArray3Variable(module, "MESH_DISPLACEMENT")
    CreateArray3Variable(module, "ROTATION")
    CreateArray3Variable(module, "VELOCITY")
    CreateArray3Variable(module, "POINT_LOAD")
    CreateArray3Variable(module, "FORCE")
    CreateArray3Variable(module, "REACTION")
    CreateArray3Variable(module, "EXTERNAL_FORCE")
    CreateArray3Variable(module, "MOMENT")
    CreateArray3Variable(module, "TORQUE")
    CreateArray3Variable(module, "NORMAL")

    CreateDoubleVariable(module, "PRESSURE")
    CreateDoubleVariable(module, "YOUNG_MODULUS")
    CreateDoubleVariable(module, "POISSON_RATIO")
    CreateDoubleVariable(module, "DOMAIN_SIZE")
    CreateDoubleVariable(module, "DENSITY")
    CreateDoubleVariable(module, "VISCOSITY")
    CreateDoubleVariable(module, "TIME")
    CreateDoubleVariable(module, "STEP")
    CreateDoubleVariable(module, "DELTA_TIME")
    CreateDoubleVariable(module, "TEMPERATURE")
    CreateDoubleVariable(module, "NODAL_MASS")
    CreateDoubleVariable(module, "NODAL_ERROR")

    CreateVectorVariable(module, "EXTERNAL_FORCES_VECTOR")
