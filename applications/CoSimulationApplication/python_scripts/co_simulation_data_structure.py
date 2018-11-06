__KRATOS_DATA_STRUCTURE__ = None

def Initialize(str_type):
    global __KRATOS_DATA_STRUCTURE__
    if(__KRATOS_DATA_STRUCTURE__ != None):
        return
    else:
        if(str_type == 'kratos'):
            __KRATOS_DATA_STRUCTURE__ = __import__('KratosMultiphysics')
        else:
            __KRATOS_DATA_STRUCTURE__ = __import__('base_co_simulation_class.co_simulation_native_data_structure')