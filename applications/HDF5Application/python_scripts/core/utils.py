'''HDF5 core utils.

license: HDF5Application/license.txt
'''


__all__ = ['DefaultSetter', 'ParametersWrapper']


from collections import Mapping
import KratosMultiphysics


class DefaultSetter(object):
    '''A utility for default settings in json.'''

    def __init__(self, settings):
        self.settings = settings

    def AddString(self, key, value):
        if not self.settings.Has(key):
            self.settings.AddEmptyValue(key).SetString(value)
        return self

    def AddInt(self, key, value):
        if not self.settings.Has(key):
            self.settings.AddEmptyValue(key).SetInt(value)
        return self

    def AddDouble(self, key, value):
        if not self.settings.Has(key):
            self.settings.AddEmptyValue(key).SetDouble(value)
        return self

    def AddBool(self, key, value):
        if not self.settings.Has(key):
            self.settings.AddEmptyValue(key).SetBool(value)
        return self

    def Add(self, key, value=KratosMultiphysics.Parameters()):
        if not self.settings.Has(key):
            if isinstance(value, str):
                self.AddString(key, value)
            elif isinstance(value, int):
                self.AddInt(key, value)
            elif isinstance(value, float):
                self.AddDouble(key, value)
            elif isinstance(value, bool):
                self.AddBool(key, value)
            else:
                self.settings.AddValue(key, value)
        return self

    def AddArray(self, key, values):
        if not self.settings.Has(key):
            self.settings.AddEmptyArray(key)
            array = self.settings[key]
            for value in values:
                array.Append(value)


class ParametersWrapper(Mapping):
    '''A pythonic wrapper to KratosMultiphysics.Parameters.

    The idea is to reduce boilerplate and improve call-site readability by
    making Parameters more pythonic without breaking existing code.
    '''

    @property
    def parameters(self):
        return self._kratos_parameters

    def __init__(self, params):
        if isinstance(params, self.__class__):
            self._kratos_parameters = params.parameters
        elif isinstance(params, str):
            self._kratos_parameters = KratosMultiphysics.Parameters(params)
        else:
            self._kratos_parameters = params

    def __getitem__(self, key):
        try:
            param = self.parameters[key]
        except:
            raise KeyError
        if param.IsString():
            value = param.GetString()
        elif param.IsDouble():
            value = param.GetDouble()
        elif param.IsBool():
            value = param.GetBool()
        elif param.IsInt():
            value = param.GetInt()
        else:
            value = self.__class__(param)
        return value

    def convert_list_to_parameters(self, list_):
        dummy_params = KratosMultiphysics.Parameters()
        dummy_params.AddEmptyList('list')
        params = dummy_params['list']
        for v in list_:
            if isinstance(v, list):
                params.Append(self.convert_list_to_parameters(v))
            else:
                params.Append(self.convert_value_to_parameters(v))
        return params

    def convert_value_to_parameters(self, value):
        params = KratosMultiphysics.Parameters()
        if isinstance(value, str):
            params.SetString(value)
        elif isinstance(value, float):
            params.SetDouble(value)
        elif isinstance(value, bool):
            params.SetBool(value)
        elif isinstance(value, int):
            params.SetInt(value)
        elif isinstance(value, KratosMultiphysics.Parameters):
            params = value
        elif isinstance(value, self.__class__):
            params = value.parameters
        else:
            raise TypeError()
        return params

    def __setitem__(self, key, value):
        if isinstance(value, list):
            params_object = self.convert_list_to_parameters(value)
        else:
            params_object = self.convert_value_to_parameters(value)
        if self.parameters.IsArray() or self.parameters.Has(key):
            self.parameters[key] = params_object
        else:
            self.parameters.AddValue(key, params_object)

    def __iter__(self):
        if self.parameters.IsArray():
            yield from range(len(self))
        else:
            yield from self.parameters.keys()

    def __len__(self):
        return self.parameters.size()

    def __str__(self):
        return self.parameters.PrettyPrintJsonString()
