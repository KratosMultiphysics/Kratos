__KRATOS_DATA_STRUCTRURE__ = None

def Initialize(str_type):
    global __KRATOS_DATA_STRUCTRURE__
    if(__KRATOS_DATA_STRUCTRURE__ != None):
        return
    else:
        if(str_type == 'kratos'):
            __KRATOS_DATA_STRUCTRURE__ = __import__('KratosMultiphysics')
        else:
            __KRATOS_DATA_STRUCTRURE__ = __import__('base_co_simulation_class.co_simulation_native_data_structure')