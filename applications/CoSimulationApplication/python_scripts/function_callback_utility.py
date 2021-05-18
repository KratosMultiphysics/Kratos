# executing a function specified as a string, similar to "kratos/utilities/python_function_callback_utility.h"

from math import *

# list of safe methods that can be used
safe_methods_list = ['acos', 'asin', 'atan', 'atan2', 'ceil', 'cos',
             'cosh', 'degrees', 'e', 'exp', 'fabs', 'floor',
             'fmod', 'frexp', 'hypot', 'ldexp', 'log', 'log10',
             'modf', 'pi', 'pow', 'radians', 'sin', 'sinh', 'sqrt',
             'tan', 'tanh']

safe_dict = {}
for k in safe_methods_list:
    safe_dict[k] = locals()[k]

def GetWordsFromString(string):
    words = [""]
    new_word = False
    for letter in string:
        if letter.isalpha():
            if new_word:
                words.append("")
                new_word = False
            words[-1]+=letter
        else:
            new_word = (words[-1]!="")

    return words

def GenericCallFunction(func_string, scope_vars=None, check=True):
    if scope_vars is not None:
        safe_dict.update(scope_vars)

    if check: # can be disabled to improve performance
        splitted_func_args = GetWordsFromString(func_string)

        for func_arg in splitted_func_args:
            if not func_arg == "" and not func_arg in safe_dict:
                raise Exception('Argument "{}" in function string "{}" was not recognized!\nOnly the following expressions can be used:\n\t{}'.format(func_arg, func_string, "\n\t".join(sorted(safe_dict.keys()))))

    return eval(func_string, {"__builtins__" : None}, safe_dict)
