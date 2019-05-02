import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyThermalInterfaceProcess(Model, settings["Parameters"])


class ApplyThermalInterfaceProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Compare against default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name" : ""
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty thermal interface model part name string. Set a valid model part name.")

        # Set the INLET flag in the thermal interface nodes and conditions
        interface_model_part = Model[settings["model_part_name"].GetString()]
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, interface_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, interface_model_part.Conditions)

        # Search for the maximum property id in the current model
        # Note that this operation is required to be performed after
        # the properties have been already set by reading the materials.json
        max_prop_id = -1
        model_part_names = Model.GetModelPartNames()
        model_part = Model.GetModelPart(settings["model_part_name"].GetString())
        for model_part_name in model_part_names:
            for prop in Model.GetModelPart(model_part_name).Properties:
                if prop.Id > max_prop_id:
                    max_prop_id = prop.Id
        model_part.GetCommunicator().MaxAll(max_prop_id)

        # Create a new property with the user defined interface parameters
        thermal_interface_prop = KratosMultiphysics.Properties(max_prop_id + 1)
        thermal_interface_prop.SetValue(KratosMultiphysics.EMISSIVITY, 0.0)
        thermal_interface_prop.SetValue(KratosMultiphysics.CONVECTION_COEFFICIENT, 0.0)

        # Set the property in the interface model part
        model_part.AddProperties(thermal_interface_prop)
        for condition in model_part.Conditions:
            condition.Properties = thermal_interface_prop
