from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
from copy import deepcopy

import KratosMultiphysics as KM

class Variable(object):
    def __init__(self, var_name, var_type, zero_val):
        self.__name = var_name
        self.__type = var_type
        self.__zero = zero_val

    def Name(self):
        return self.__name

    def Type(self):
        return self.__type

    def Zero(self):
        # copy to make sure that nothing is referenced wrong
        return deepcopy(self.__zero)

    def __hash__(self):
        return hash(self.__name)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return 'Variable "{}" of type "{}"'.format(self.__name, self.__type)

class VariableComponent(Variable):
    def __init__(self, var_name, source_variable, component_index):
        super(VariableComponent, self).__init__(var_name, "Component", 0.0)
        self.__source_variable = source_variable
        self.__component_index = component_index

    def GetSourceVariable(self):
        return self.__source_variable

    def GetComponentIndex(self):
        return self.__component_index

    def __str__(self):
        return 'Variable-Component "{}"'.format(self.Name())


def CreateDoubleVariable(module, name):
    if name in KM._registered_variables:
        raise NameError('Variable "{}" exists already!'.format(name))

    KM._registered_variables[name] = Variable(name, "Double", 0.0)
    setattr(module, name, KM._registered_variables[name])

def CreateComponentVariable(module, name, source_variable, component_index):
    if name in KM._registered_variables:
        raise NameError('Variable "{}" exists already!'.format(name))
    KM._registered_variables[name] = VariableComponent(name, source_variable, component_index)
    setattr(module, name, KM._registered_variables[name])

def CreateArray3Variable(module, name):
    if name in KM._registered_variables:
        raise NameError('Variable "{}" exists already!'.format(name))

    array_var = Variable(name, "Array", [0.0, 0.0, 0.0])
    KM._registered_variables[name] = array_var
    setattr(module, name, KM._registered_variables[name])

    for i_comp, comp in enumerate(["X", "Y", "Z"]):
        CreateComponentVariable(module, name + "_" + comp, array_var, i_comp)

def CreateVectorVariable(module, name):
    if name in KM._registered_variables:
        raise NameError('Variable "{}" exists already!'.format(name))
    KM._registered_variables[name] = Variable(name, "Vector", [])
    setattr(module, name, KM._registered_variables[name])
