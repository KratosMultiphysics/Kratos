import KratosMultiphysics as KM

class DeprecationManager:

    @classmethod
    def HasDeprecatedVariable(cls, context_string, settings, old_variable_name, new_variable_name):

        if settings.Has(old_variable_name):
            if not settings.Has(new_variable_name):
                KM.Logger.PrintWarning(context_string,
                                       '\n\x1b[1;31m[DEPRECATED INPUT PARAMETERS] \x1b[0m' + '\''
                                       + old_variable_name + '\' is deprecated; use \''
                                       + new_variable_name + '\' instead.')
                return True
            else:
                raise NameError('Conflicting input variable names: Both the deprecated variable \''
                                + old_variable_name + '\' and its current standard replacement \''
                                + new_variable_name + '\' were found. Please, remove \''
                                + old_variable_name + '\'.')
        return False

    @classmethod
    def ReplaceDeprecatedVariableName(cls, settings, old_variable_name, new_variable_name):
        settings.AddEmptyValue(new_variable_name)

        if settings[old_variable_name].IsInt():
            settings[new_variable_name].SetInt(settings[old_variable_name].GetInt())
        elif settings[old_variable_name].IsDouble():
            settings[new_variable_name].SetDouble(settings[old_variable_name].GetDouble())
        elif settings[old_variable_name].IsString():
            settings[new_variable_name].SetString(settings[old_variable_name].GetString())
        elif settings[old_variable_name].IsBool():
            settings[new_variable_name].SetBool(settings[old_variable_name].GetBool())
        else:
            pass

        settings.RemoveValue(old_variable_name)