from __future__ import print_function, absolute_import, division

from utils.constants import ctab


class VariableCreator(object):
    def __init__(self, name, vtype, is3D=False):
        ''' Creates a variable for an application

            Input
            -----
            - name: string
                name of the variable

            - vtype: string
                type of the variable

            - is3D:: boolean
                determines if the variable is vectorial(True) or scalar(False, default)
                NOTE: This will be could be replaced by VariableCreator3D at some point.

        '''

        self.defineString = 'KRATOS_DEFINE_VARIABLE( {type}, {name} )\n'.format(type=vtype, name=name)
        self.createString = 'KRATOS_CREATE_VARIABLE( {type}, {name} )\n'.format(type=vtype, name=name)
        self.registerString = ctab + 'KRATOS_REGISTER_VARIABLE( {name} )\n'.format(name=name)
        self.registerPythonString = ctab + 'KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, {name} )\n'.format(name=name)

        # String changes if is a 3D variable
        if is3D:
            self.defineString = 'KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
            self.createString = 'KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
            self.registerString = ctab + 'KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
            self.registerPythonString = ctab + 'KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, {name} )\n'.format(name=name)


class VariableCreator3D(object):
    def __init__(self, name):
        ''' Creates a 3D variable for an application.
            All 3D variables are "double" by definition

            Input
            -----
            - name: string
                name of the variable
        '''

        self.defineString = 'KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
        self.createString = 'KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
        self.registerString = ctab + 'KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
        self.registerPythonString = ctab + 'KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
