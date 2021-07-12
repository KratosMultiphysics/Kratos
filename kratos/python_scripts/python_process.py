import KratosMultiphysics

class PythonProcess(KratosMultiphysics.Process):
    """This process is a dummy python process which derives directly from the C++ Process class, so it has all the methods relative to the the Process clas
    """
    
    def __init__(self):
        """ The default constructor of the class
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        KratosMultiphysics.Process.__init__(self)
        
def Factory(settings, Model):
    return PythonProcess()
