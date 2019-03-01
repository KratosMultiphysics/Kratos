import json


## To pass the parameters to different objects.
class Parameters(object):
    def __init__(self, text_stream):
        if(isinstance(text_stream, str)):
            self.parameters = json.loads(text_stream)
        elif(isinstance(text_stream, dict) or isinstance(text_stream, list) ):
            self.parameters = text_stream
        #super(Parameters, self).__init__(self.parameters)
        self.Initialize()
        self.count = -1

    def Initialize(self):
        pass
        #
        # if(not isinstance(self.parameters, list)):
            #super(Parameters, self).__init__(self.parameters)

    def size(self):
        return len(self.parameters)

    def ValidateAndAssignDefaults(self, defaults, recursive=False):
        for key, val in self.parameters.items():
            # check if the current entry also exists in the defaults
            if not key in defaults.keys():
                err_msg  = 'The item with name "' + key + '" is present in this '
                err_msg += 'self.parameters\nbut NOT in the defaults!\n'
                err_msg += 'self.parameters are:\n'
                err_msg += json.dumps(self.parameters, indent=4)
                err_msg += '\ndefaults are:\n'
                err_msg += json.dumps(defaults, indent=4)
                raise Exception(err_msg)

            # check if the type is the same in the defaults
            if type(self.parameters[key]) != defaults[key].Type():
                err_msg  = 'The type of the item with name "' + key + '" (type: "'
                err_msg += str(type(self.parameters[key]).__name__)+'") in this '
                err_msg += 'self.parameters\nis NOT the same as in the defaults (type: "'
                err_msg += str(defaults[key].Type().__name__)+'")!\n'
                err_msg += 'self.parameters are:\n'
                err_msg += json.dumps(self.parameters, indent=4)
                err_msg += '\ndefaults are:\n'
                err_msg += json.dumps(defaults, indent=4)
                raise Exception(err_msg)

        # loop the defaults and add the missing entries
        for key_d, val_d in defaults.items():
            if key_d not in self.parameters: # add the default in case the setting is not present
                self.parameters[key_d] = val_d
            elif recursive and type(val_d) is dict:
                RecursivelyValidateAndAssignDefaults(val_d, self.parameters[key_d])

        self.Initialize()

    def RecursivelyValidateAndAssignDefaults(self, defaults):
        ValidateAndAssignDefaults(defaults, self.parameters, recursive=True)

    def __getitem__(self, key):
        a = Parameters("{}")
        if( isinstance(self.parameters[key], dict) or isinstance(self.parameters[key], list)):
            a.parameters = self.parameters[key]
            a.Initialize()
            return a
        else:
            a.parameters = {"a_py_kratos":self.parameters[key]}
            a.Initialize()
            #return self.parameters[key]
            return a

    def __iter__(self):
        self.count = -1
        return self

    def __next__(self):
        self.count = self.count + 1
        if(self.count<len(self.parameters)):
            return Parameters(self.parameters[self.count])
        else:
            raise StopIteration

    def __contains__(self, item):
        return (item in self.parameters)

    def keys(self):
        if(isinstance(self.parameters, dict)):
            return self.parameters.keys()

    def items(self):
        if(isinstance(self.parameters, dict)):
            #key, value = self.parameters.items()
            return self.parameters.items()
            #return (key, Parameters(value))

    def GetString(self):
        return str(self.parameters['a_py_kratos'])
    def GetDouble(self):
        return float(self.parameters['a_py_kratos'])
    def GetInt(self):
        return int(self.parameters['a_py_kratos'])
    def IsInt(self):
        return isinstance(self.parameters['a_py_kratos'], int)
    def IsDouble(self):
        return isinstance(self.parameters['a_py_kratos'], float)
    def IsString(self):
        return isinstance(self.parameters['a_py_kratos'], str)
    def Has(self, param):
        return (param in self.parameters)

    def AddValue(self, key, value):
        if(isinstance(self.parameters, dict)):
            if(key not in self.parameters):
                self.parameters[key] = value
            else:
                RuntimeError("Key already exists")

    def Type(self):
        if( isinstance(self.parameters, list) ):
            return type(self.parameters)
        elif(isinstance(self.parameters, dict)) :
            if('a_py_kratos' in self.parameters):
                return type(self.parameters['a_py_kratos'])
            else:
                return type(self.parameters)
