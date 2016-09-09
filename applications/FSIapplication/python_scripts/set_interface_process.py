from KratosMultiphysics import * 

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetInterfaceProcess(Model, settings["Parameters"])

class SetInterfaceProcess(Process):
    def __init__(self, Model, settings):
        
        Process.__init__(self)
        
        interface_model_part = Model[settings["model_part_name"].GetString()]
        
        if settings["variable_name"].GetString() == "STRUCTURE_INTERFACE":
            for node in interface_model_part.Nodes:
                # Set the INTERFACE flag
                node.Set(INTERFACE, True)
        
        elif settings["variable_name"].GetString() == "FLUID_INTERFACE":
            zero_vect = [0,0,0]
            
            for node in interface_model_part.Nodes:
                # Set the INTERFACE flag
                node.Set(INTERFACE, True)
                # Fix the DISPLACEMENT and initialize it to zero at the interface
                node.Fix(DISPLACEMENT_X)
                node.Fix(DISPLACEMENT_Y)
                node.Fix(DISPLACEMENT_Z)
                node.SetSolutionStepValue(DISPLACEMENT,0,zero_vect)
                
        
