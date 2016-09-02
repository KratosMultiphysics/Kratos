from KratosMultiphysics import * 

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetInterfaceProcess(Model, settings["Parameters"])

class SetInterfaceProcess(Process):
    def __init__(self, Model, settings):
        
        Process.__init__(self)
        
        interface_model_part = Model[settings["model_part_name"].GetString()]
        
        for node in interface_model_part.Nodes:
            node.Set(INTERFACE, True)
