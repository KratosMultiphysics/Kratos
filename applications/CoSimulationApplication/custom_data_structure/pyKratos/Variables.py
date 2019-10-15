from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
from copy import deepcopy

class Variable(object):
    """
    Variables are global objects which are a part of the
    KratosMultiphysics module.

    They can be accessed in different ways, see code
    example below:

        import KratosMultiphysics as KM
        var = vars(KM)[disp]
        var = KM.__dict__[disp]
        var = KM.KratosGlobals.GetVariable(disp)

    Variables must be created in the module itself,
    i.e. in the Variables.py file. As far as I know,
    no new Variables can be created on-the-go.
    Therefore, if you need a variable which has not
    been defined yet, add it to this file.

    For 1D variables, CreateDoubleVariable is used,
    for 3D vector variables, CreateArray3Variable is used.
    Not sure what CreateVectorVariable is for...
    """

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


# TODO maybe separate the class-declarations from the variable declarations?
def CreateDoubleVariable(name):
    if name in globals():
        raise NameError('Variable "{}" exists already!'.format(name))
    globals()[name] = Variable(name, "Double", 0.0)

def CreateComponentVariable(name, source_variable, component_index):
    if name in globals():
        raise NameError('Variable "{}" exists already!'.format(name))
    globals()[name] = VariableComponent(name, source_variable, component_index)

def CreateArray3Variable(name):
    if name in globals():
        raise NameError('Variable "{}" exists already!'.format(name))

    array_var = Variable(name, "Array", [0.0, 0.0, 0.0])
    globals()[name] = array_var

    for i_comp, comp in enumerate(["X", "Y", "Z"]):
        CreateComponentVariable(name + "_" + comp, array_var, i_comp)

def CreateVectorVariable(name):
    if name in globals():
        raise NameError('Variable "{}" exists already!'.format(name))
    globals()[name] = Variable(name, "Vector", [])

CreateArray3Variable("DISPLACEMENT")
CreateArray3Variable("MESH_DISPLACEMENT")
CreateArray3Variable("ROTATION")
CreateArray3Variable("VELOCITY")
CreateArray3Variable("POINT_LOAD")
CreateArray3Variable("FORCE")
CreateArray3Variable("REACTION")
CreateArray3Variable("EXTERNAL_FORCE")
CreateArray3Variable("TORQUE")
CreateArray3Variable("NORMAL")

CreateDoubleVariable("PRESSURE")
CreateDoubleVariable("YOUNG_MODULUS")
CreateDoubleVariable("POISSON_RATIO")
CreateDoubleVariable("DOMAIN_SIZE")
CreateDoubleVariable("DENSITY")
CreateDoubleVariable("VISCOSITY")
CreateDoubleVariable("TIME")
CreateDoubleVariable("DELTA_TIME")
CreateDoubleVariable("TEMPERATURE")
CreateDoubleVariable("NODAL_MASS")
CreateDoubleVariable("NODAL_ERROR")

CreateVectorVariable("EXTERNAL_FORCES_VECTOR")

CreateDoubleVariable("AREA")
