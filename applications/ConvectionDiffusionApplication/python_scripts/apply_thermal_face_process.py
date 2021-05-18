import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyThermalFaceProcess(Model, settings["Parameters"])


class ApplyThermalFaceProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Compare with default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "ambient_temperature": 0.0,
            "add_ambient_radiation": false,
            "emissivity": 0.0,
            "add_ambient_convection": false,
            "convection_coefficient": 0.0,
            "interval": [0.0,1e30]
        }''')

        settings.ValidateAndAssignDefaults(default_settings)

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty thermal face model part name string. Set a valid model part name.")
        if not settings["ambient_temperature"].IsDouble():
            raise Exception("Ambient temperature must be a constant value. Time dependency not implemented yet.")

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
        max_prop_id = model_part.GetCommunicator().GetDataCommunicator().MaxAll(max_prop_id)

        # Create a new property with the user defined interface parameters
        thermal_interface_prop = KratosMultiphysics.Properties(max_prop_id + 1)
        emissivity = settings["emissivity"].GetDouble() if settings["add_ambient_radiation"].GetBool() else 0.0
        convection_coefficient = settings["convection_coefficient"].GetDouble() if settings["add_ambient_convection"].GetBool() else 0.0
        ambient_temperature = settings["ambient_temperature"].GetDouble() if (settings["add_ambient_convection"].GetBool() or settings["add_ambient_radiation"]) else 0.0
        thermal_interface_prop.SetValue(KratosMultiphysics.EMISSIVITY, emissivity)
        thermal_interface_prop.SetValue(KratosMultiphysics.CONVECTION_COEFFICIENT, convection_coefficient)
        thermal_interface_prop.SetValue(KratosMultiphysics.AMBIENT_TEMPERATURE, ambient_temperature)

        # Set the new property in the thermal face model part
        model_part.AddProperties(thermal_interface_prop)
        for condition in model_part.Conditions:
            condition.Properties = thermal_interface_prop