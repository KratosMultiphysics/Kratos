'''HDF5 core utils.

license: HDF5Application/license.txt
'''
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
