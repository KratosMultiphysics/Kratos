import KratosMultiphysics

class PythonOperation(object):
    def __init__(self):
        KratosMultiphysics.Logger.PrintWarning("PythonOperation", "This is a fake class to be used while developing. Must be implemented in C++.")

    def Execute(self):
        KratosMultiphysics.Logger.PrintWarning("PythonOperation", "Calling the fake PythonOperation Execute().")

def Factory(settings, model):
    return PythonOperation()
