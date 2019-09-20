from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

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

def GenericCallFunction(func_string, scope_vars):
    safe_dict.update(scope_vars)
    return eval(func_string, {"__builtins__" : None}, safe_dict)
