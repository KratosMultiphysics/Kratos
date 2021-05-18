'''HDF5 core utils.

license: HDF5Application/license.txt
'''


__all__ = ['ParametersWrapper']


from collections.abc import Mapping


import KratosMultiphysics


class ParametersWrapper(Mapping):
    '''A pythonic wrapper to KratosMultiphysics.Parameters.

    The idea is to reduce boilerplate and improve call-site readability by
    making Parameters more pythonic without breaking existing code.
    '''

    def __init__(self, params="{}"):
        if isinstance(params, self.__class__):
            self._parameters = params._parameters
        elif isinstance(params, str):
            self._parameters = KratosMultiphysics.Parameters(params)
        else:
            self._parameters = params

    def __getitem__(self, key):
        try:
            param = self._parameters[key]
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

    def _convert_list_to_parameters(self, list_):
        dummy_params = KratosMultiphysics.Parameters()
        dummy_params.AddEmptyList('list')
        params = dummy_params['list']
        for v in list_:
            if isinstance(v, list):
                params.Append(self._convert_list_to_parameters(v))
            else:
                params.Append(self._convert_value_to_parameters(v))
        return params

    def _convert_value_to_parameters(self, value):
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
            params = value._parameters
        else:
            raise TypeError()
        return params

    def __setitem__(self, key, value):
        if isinstance(value, list):
            params_object = self._convert_list_to_parameters(value)
        else:
            params_object = self._convert_value_to_parameters(value)
        if self._parameters.IsArray() or self._parameters.Has(key):
            self._parameters[key] = params_object
        else:
            self._parameters.AddValue(key, params_object)

    def __iter__(self):
        '''Return an iterator over keys or indices.

        This iterates over keys or indices, depending on if the Parameters
        instance is an array.
        '''
        if self._parameters.IsArray():
            yield from range(len(self))
        else:
            yield from self._parameters.keys()

    def __len__(self):
        return self._parameters.size()

    def __str__(self):
        return self._parameters.PrettyPrintJsonString()

    def __getattr__(self, name):
        '''Expose methods defined in KratosMultiphysics.Parameters.'''
        return getattr(self._parameters, name)

    def Get(self):
        return self._parameters

    def SetDefault(self, key, value=KratosMultiphysics.Parameters()):
        if key not in self:
            self[key] = value
