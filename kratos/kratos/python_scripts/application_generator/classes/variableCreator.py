from utils.constants import ctab


class VariableCreator(object):
    def __init__(self, name, vtype, is3D=False):

        self.defineString = 'KRATOS_DEFINE( {type}, {name} )\n'.format(type=vtype, name=name)
        self.createString = 'KRATOS_CREATE( {type}, {name} )\n'.format(type=vtype, name=name)
        self.registerString = ctab + 'KRATOS_REGISTER( {name} )\n'.format(name=name)
        self.registerPythonString = ctab + 'KRATOS_REGISTER_IN_PYTHON_VARIABLE( {name} )\n'.format(name=name)

        # String changes if is a 3D variable
        if is3D:
            self.defineString = 'KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
            self.createString = 'KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
            self.registerString = ctab + 'KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
            self.registerPythonString = ctab + 'KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( {name} )\n'.format(name=name)
