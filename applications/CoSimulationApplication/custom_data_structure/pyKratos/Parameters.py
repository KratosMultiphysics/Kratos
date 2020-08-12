import json

class Parameters(object):

    def __init__(self, input_stream):
        if isinstance(input_stream, str):
            # Constructing from a stream
            self.param = json.loads(input_stream)
        else:
            raise TypeError('"Parameters" can only be constructed from a string in json format')

    def size(self):
        self.__CheckIfArray("size")
        return len(self.param)

    def PrettyPrintJsonString(self):
        return json.dumps(self.param, indent=4)

    def ValidateAndAssignDefaults(self, defaults, recursive=False):
        self.__CheckIfParameter(defaults)

        for key in self.keys():
            # check if the current entry also exists in the defaults
            if not defaults.Has(key):
                err_msg  = 'The item with name "' + key + '" is present in this '
                err_msg += 'settings\nbut NOT in the defaults!\n'
                err_msg += 'settings are:\n'
                err_msg += self.PrettyPrintJsonString()
                err_msg += '\ndefaults are:\n'
                err_msg += defaults.PrettyPrintJsonString()
                raise Exception(err_msg)

            # check if the type is the same in the defaults
            if type(self.param[key]) != type(defaults[key].param):
                err_msg  = 'The type of the item with name "' + key + '" (type: "'
                err_msg += str(type(self.param[key]).__name__)+'") in this '
                err_msg += 'settings\nis NOT the same as in the defaults (type: "'
                err_msg += str(type(defaults[key].param).__name__)+'")!\n'
                err_msg += 'settings are:\n'
                err_msg += self.PrettyPrintJsonString()
                err_msg += '\ndefaults are:\n'
                err_msg += defaults.PrettyPrintJsonString()
                raise Exception(err_msg)

        # loop the defaults and add the missing entries
        for key_d, val_d in defaults.items():
            if not self.Has(key_d): # add the default in case the setting is not present
                self.AddValue(key_d, val_d)
            elif recursive and val_d.IsSubParameter():
                self.__CreateParameters(self.param[key_d]).RecursivelyValidateAndAssignDefaults(val_d)

    def RecursivelyValidateAndAssignDefaults(self, defaults):
        self.ValidateAndAssignDefaults(defaults, recursive=True)

    def AddMissingParameters(self, missing_param):
        self.__CheckIfParameter(missing_param)

        for key, val in missing_param.items():
            if not self.Has(key):
                self.AddValue(key, val)

    def Has(self, key):
        return key in self.param

    def __getitem__(self, key):
        if self.IsSubParameter():
            if not self.Has(key):
                raise KeyError('Error: Getting a value that does not exist. entry string :', key)
            return self.__CreateParameters(self.param[key])
        elif self.IsArray():
            if key >= self.size():
                raise IndexError("Index exceeds array size ({}). Index value is : {}".format(self.size(), key))
            return self.__CreateParameters(self.param[key])
        else:
            raise TypeError("This object cannot be accessed by index")

    def __repr__(self):
        return "Python Parameters Object " + self.PrettyPrintJsonString()

    def __CreateParameters(self, content):
        new_param = Parameters("{}")
        new_param.param = content
        return new_param

    def keys(self):
        self.__CheckIfSubParameter("keys")
        return list(self.param.keys())

    def values(self):
        self.__CheckIfSubParameter("values")
        return [self.__CreateParameters(p) for p in self.param.values()]

    def items(self):
        self.__CheckIfSubParameter("items")
        return [tup for tup in zip(self.keys(), self.values())]

    def AddValue(self, key, other_param):
        self.__CheckIfSubParameter("AddValue")
        self.__CheckIfParameter(other_param)
        if other_param == self:
            raise Exception("The object cannot be added to itself!")
        if self.Has(key):
            print('Warning, key "{}" exists already and will be overwritten!'.format(key))

        self.param[key] = other_param.param

    def AddEmptyValue(self, key):
        self.__CheckIfSubParameter("AddEmptyValue")
        if not self.Has(key):
            self.param[key] = self.__CreateParameters(None)
        return self.param[key]

    def RemoveValue(self, key):
        self.__CheckIfSubParameter("RemoveValue")
        if self.Has(key): # removing a non-existing key does not throw in Kratos
            self.param.pop(key)


    def __CheckIfSubParameter(self, fct_name):
        if not self.IsSubParameter():
            raise TypeError('"{}" can only be used if the value is of Parameter type'.format(fct_name))

    def __CheckIfArray(self, fct_name):
        if not self.IsArray():
            raise TypeError('"{}" can only be used if the value is of Array type'.format(fct_name))

    def __CheckIfParameter(self, obj):
        if not type (obj) == Parameters:
            raise TypeError('Input is not of type "Parameters", but of type "{}" from module "{}"!'.format(obj.__class__.__name__, obj.__module__))


    #########################
    ##### IsXXX Methods #####
    #########################

    def IsInt(self):
        return self.__Is(int)

    def IsDouble(self):
        return self.__Is(float)

    def IsNumber(self):
        return self.IsInt() or self.IsDouble()

    def IsBool(self):
        return self.__Is(bool)

    def IsString(self):
        return self.__Is(str)

    def IsArray(self):
        return self.__Is(list)

    def IsSubParameter(self):
        return self.__Is(dict)

    def IsNull(self):
        return self.__Is(type(None))


    def __Is(self, exp_type):
        return type(self.param) == exp_type


    ##########################
    ##### GetXXX Methods #####
    ##########################

    def GetInt(self):
        return self.__Get(self.IsNumber, "number")

    def GetDouble(self):
        return self.__Get(self.IsNumber, "number")

    def GetBool(self):
        return self.__Get(self.IsBool, "bool")

    def GetString(self):
        return self.__Get(self.IsString, "string")

    def GetStringArray(self):
        self.__CheckIfArray("GetStringArray")
        return [self[i].GetString() for i in range(self.size())]


    def __Get(self, cmp_fct, exp_type_str):
        if not cmp_fct():
            raise TypeError("Argument must be a {}!".format(exp_type_str))
        return self.param



    ##########################
    ##### SetXXX Methods #####
    ##########################

    def SetInt(self, val):
        self.__Set(val)

    def SetDouble(self):
        self.__Set(val)

    def SetBool(self):
        self._Set(val)

    def SetString(self):
        self.__Set(val)


    def __Set(self, val):
        self.param = val
