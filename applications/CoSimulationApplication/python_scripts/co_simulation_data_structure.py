__DATA_STRUCTURE__ = None

def Initialize(str_type):
    global __DATA_STRUCTURE__
    if(__DATA_STRUCTURE__ != None):
        return
    else:
        if(str_type == 'kratos'):
            __DATA_STRUCTURE__ = __import__('KratosMultiphysics')
        elif(str_type == 'pyKratos'):
            __DATA_STRUCTURE__ = __import__('pyKratos')