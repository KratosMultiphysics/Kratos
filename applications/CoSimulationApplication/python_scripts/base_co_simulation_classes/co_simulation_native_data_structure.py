import json


## To pass the parameters to different objects.
class Parameters(dict):
    __init__(self, text_stream):
        self.parameters = json.load(test_stream)
        super(Parameters, self).__init__(self.parameters)

    def size(self):
        return len(self.parameters)

    # TODO: make it recursive by default
    def ValidateAndAssignDefaults(self, default, warnUnused=True):
        output = dict()
        for key in default.keys():
            if key in self.parameters.keys():
                value = self.parameters[key]
                if type(default[key]) == type(value) or default[key] == type(value):
                    output[key] = value
                else:
                    error =  "ERROR: Key ", key, "has not the same type as the default value! "
                    input_value = "INPUT: ", type(value)
                    default_value = "DEFAULT: " + str(type(default[key])) + str(type(value))
                    raise ValueError(error, input_value, default_value)

            else:
                if type(default[key]) == type:
                    error =  "ERROR: mandatory setting ", key, ":", default[key], "is missing! "
                    raise ValueError(error)
                else:
                    output[key] = default[key]

            for inputSetting in self.parameters.keys():
                if inputSetting not in default.keys():
                    value = self.parameters[inputSetting]
                    if type(value) == str:
                        value = value.encode('utf-8')
                    output[inputSetting] = value

                    if warnUnused:
                        warning_msg = "UNUSED setting in self.parameters: " + str(inputSetting) + ":" + str(value)
                        warnings.warn(warning_msg, Warning)

        return output

## To store all the meshes imported or to be exported.
class Model(object):
    pass

# ModelPart class
class ModelPart(object):
    pass