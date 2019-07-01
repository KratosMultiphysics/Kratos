from math import *

# list of safe methods that can be used
safe_list = ['acos', 'asin', 'atan', 'atan2', 'ceil', 'cos',
             'cosh', 'degrees', 'e', 'exp', 'fabs', 'floor',
             'fmod', 'frexp', 'hypot', 'ldexp', 'log', 'log10',
             'modf', 'pi', 'pow', 'radians', 'sin', 'sinh', 'sqrt',
             'tan', 'tanh']


safe_dict = {}
for k in safe_list:
    safe_dict[k] = locals()[k]

def GenericCallFunction(func_string, scope_vars):
    safe_dict.update(scope_vars)
    return eval(func_string, {"__builtins__" : None}, safe_dict)
