from KratosMultiphysics import *

class KratosProcessFactory(object):
    def __init__(self, Model):
        self.Model = Model #model is a place
        
    def ConstructListOfProcesses( self, process_list ):
        constructed_processes = []
        for item in process_list:
            if( "implemented_in_python" in item): #check if implemented in python or in c++
                if item["implemented_in_python"] == True: #here we treat the case of python implemented processes
                    module = __import__(item["implemented_in_module"])
                    interface_file = __import__(item["implemented_in_file"])
                    p = interface_file.Factory(item, self.Model)
                    constructed_processes.append( p )
                    
                else: #here we create c++ processes
                    module = __import__(item["implemented_in_module"])
                    interface_file = __import__(item["implemented_in_file"])
                    p = interface_file.Factory(item, self.Model)
                    constructed_processes.append( p )
                    
            else:
                raise Exception( "process defined by settings ", item, " does not define if implemented in python or in c++")
            
    
        return constructed_processes