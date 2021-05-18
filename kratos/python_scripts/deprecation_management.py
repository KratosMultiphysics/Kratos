import KratosMultiphysics as KM

class DeprecationManager:
    """
    This class is intended to encapsulate common operations that may be needed
    when dealing with deprecated input variable names. Its original purpose is
    the management of json-type input, although it may be extended to other input
    types.

    The basic goals that inspired this encapsulation are:
        1. Avoid repeating too much code.
        2. Make these operations more robust.
        3. Facilitate the evolution of terminology by concentrating the changes to be made
        around a single object.
        4. Encourage better legacy-terminology handling, including appropriate error/warning handling.
        5. Providing standard tools useful in the creation of conventional protocols for the modification
        of the API.
    """
    @staticmethod
    def HasDeprecatedVariable(context_string, parameters, old_variable_name, new_variable_name):
        """ Check if a given deprecated variable is present in a given parameters object (Parameters object)"""
        if parameters.Has(old_variable_name):
            if not parameters.Has(new_variable_name):
                KM.Logger.PrintWarning(context_string,
                                       '\n\x1b[1;31m[DEPRECATED INPUT PARAMETERS] \x1b[0m' + '\''
                                       + old_variable_name + '\' is deprecated; use \''
                                       + new_variable_name + '\' instead.')
                return True
            else: # Error: the new and deprecated variables cannot coexist
                raise NameError('Conflicting input variable names: Both the deprecated variable \''
                                + old_variable_name + '\' and its current standard replacement \''
                                + new_variable_name + '\' were found. Please, remove \''
                                + old_variable_name + '\'.')
        return False

    @staticmethod
    def ReplaceDeprecatedVariableName(parameters, old_variable_name, new_variable_name):
        """ Replace a key by another. The old key is assumed to be present."""
        parameters.AddEmptyValue(new_variable_name)

        if parameters[old_variable_name].IsInt():
            parameters[new_variable_name].SetInt(parameters[old_variable_name].GetInt())
        elif parameters[old_variable_name].IsDouble():
            parameters[new_variable_name].SetDouble(parameters[old_variable_name].GetDouble())
        elif parameters[old_variable_name].IsString():
            parameters[new_variable_name].SetString(parameters[old_variable_name].GetString())
        elif parameters[old_variable_name].IsBool():
            parameters[new_variable_name].SetBool(parameters[old_variable_name].GetBool())
        else:
            pass

        parameters.RemoveValue(old_variable_name)