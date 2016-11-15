import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetInterfaceProcess(Model, settings["Parameters"])

class SetInterfaceProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "mesh_id"         : 0,
            "model_part_name" : "CHOOSE_FLUID_OR_STRUCTURE_INTERFACE_MODELPART_NAME",
            "variable_name"   : "CHOOSE_BETWEEN_STRUCTURE_INTERFACE_OR_FLUID_INTERFACE"
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);
        
        interface_model_part = Model[settings["model_part_name"].GetString()]
        
        if settings["variable_name"].GetString() == "STRUCTURE_INTERFACE":
            for node in interface_model_part.Nodes:
                # Set the INTERFACE flag
                node.Set(KratosMultiphysics.INTERFACE, True)
        
        elif settings["variable_name"].GetString() == "FLUID_INTERFACE":
            zero_vect = [0,0,0]
            
            for node in interface_model_part.Nodes:
                # Set the INTERFACE flag
                node.Set(KratosMultiphysics.INTERFACE, True)
                # Fix the DISPLACEMENT and initialize it to zero at the interface
                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,zero_vect)
