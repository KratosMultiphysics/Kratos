import json


## To pass the parameters to different objects.
class Parameters(dict):
    __init__(self, text_stream):
        self.parameters = json.load(test_stream)
        super(Parameters, self).__init__(self.parameters)

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
            if type(self.parameters[key]) != type(defaults[key]):
                err_msg  = 'The type of the item with name "' + key + '" (type: "'
                err_msg += str(type(self.parameters[key]).__name__)+'") in this '
                err_msg += 'self.parameters\nis NOT the same as in the defaults (type: "'
                err_msg += str(type(defaults[key]).__name__)+'")!\n'
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

    def RecursivelyValidateAndAssignDefaults(self, defaults):
        ValidateAndAssignDefaults(defaults, self.parameters, recursive=True)
