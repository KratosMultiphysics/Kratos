import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AutomaticDirichletConditionsProcess(model, settings["Parameters"])

## This process fix the variables at the skin
class AutomaticDirichletConditionsProcess(KM.Process):

    def __init__(self, model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"             : "main_model_part",
                "scalar_variables_list"       : [],
                "scalar_value"                : 0.0,
                "scalar_constraint"           : true,
                "vector_variables_list"       : [],
                "vector_value"                : [0.0, 0.0, 0.0],
                "vector_constraint"           : [true, true, true],
                "auxiliary_model_part_name"   : "skin_model_part",
                "auxiliary_condition_name"    : "Condition",
                "skin_parts_to_exclude"       : [],
                "allow_outwards_flow_on_land" : true,
                "sea_water_level"             : 0.0
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)
        self.settings = settings
        self.model = model

        self._CheckIfParameterIsRelativeSubModelPartName(self.settings["auxiliary_model_part_name"])
        for item in self.settings["skin_parts_to_exclude"]:
            self._CheckIfParameterIsRelativeSubModelPartName(item)

    def ExecuteInitialize(self):
        model_part_name = self.settings["model_part_name"].GetString()
        model_part = self.model[model_part_name]

        skin_detection_settings = KM.Parameters()
        skin_detection_settings.AddValue("name_auxiliar_model_part", self.settings["auxiliary_model_part_name"])
        skin_detection_settings.AddValue("name_auxiliar_condition", self.settings["auxiliary_condition_name"])
        skin_detection_settings.AddValue("list_model_parts_to_assign_conditions", self.settings["skin_parts_to_exclude"])
        KM.SkinDetectionProcess2D(model_part, skin_detection_settings).Execute()

        skin_model_part_name = model_part_name + '.' + self.settings["auxiliary_model_part_name"].GetString()
        skin_model_part = self.model[skin_model_part_name]

        # We remove the new interface condition from the bondaries with other conditions
        for name in self.settings["skin_parts_to_exclude"].GetStringArray():
            skin_to_exclude = model_part_name + '.' + name
            self.model[skin_to_exclude].RemoveConditionsFromAllLevels(KM.INTERFACE)
            KM.VariableUtils().SetFlag(KM.INTERFACE, False, self.model[skin_to_exclude].Nodes)
        skin_model_part.RemoveNodes((KM.INTERFACE).AsFalse())

        # Here we exclude the interface which will has an outward flow
        if self.settings["allow_outwards_flow_on_land"].GetBool():
            dimension = 2
            KM.NormalCalculationUtils().CalculateOnSimplex(skin_model_part, dimension)
            KM.ComputeNonHistoricalNodalGradientProcess(model_part, SW.TOPOGRAPHY, SW.TOPOGRAPHY_GRADIENT, KM.NODAL_AREA).Execute()
            SW.ShallowWaterUtilities().IdentifySolidBoundary(skin_model_part, self.settings["sea_water_level"].GetDouble(), KM.SOLID)
            skin_model_part.RemoveNodes((KM.SOLID).AsFalse())
            skin_model_part.RemoveConditionsFromAllLevels((KM.SOLID).AsFalse())
        else:
            KM.VariableUtils().SetFlag(KM.SOLID, True, skin_model_part.Nodes)
            KM.VariableUtils().SetFlag(KM.SOLID, True, skin_model_part.Conditions)

        self.processes = []

        from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess
        for variable in self.settings["scalar_variables_list"].GetStringArray():
            scalar_settings = KM.Parameters()
            scalar_settings.AddEmptyValue("model_part_name").SetString(skin_model_part_name)
            scalar_settings.AddEmptyValue("variable_name").SetString(variable)
            scalar_settings.AddValue("value", self.settings["scalar_value"])
            scalar_settings.AddValue("constrained", self.settings["scalar_constraint"])
            self.processes.append(AssignScalarVariableProcess(self.model, scalar_settings))

        from KratosMultiphysics.assign_vector_variable_process import AssignVectorVariableProcess
        for variable in self.settings["vector_variables_list"].GetStringArray():
            vector_settings = KM.Parameters()
            vector_settings.AddEmptyValue("model_part_name").SetString(skin_model_part_name)
            vector_settings.AddEmptyValue("variable_name").SetString(variable)
            vector_settings.AddValue("value", self.settings["vector_value"])
            vector_settings.AddValue("constrained", self.settings["vector_constraint"])
            self.processes.append(AssignVectorVariableProcess(self.model, vector_settings))

    def ExecuteBeforeSolutionLoop(self):
        for process in self.processes:
            process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        for process in self.processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for process in self.processes:
            process.ExecuteFinalizeSolutionStep()

    @staticmethod
    def _CheckIfParameterIsRelativeSubModelPartName(param):
        name = param.GetString()
        search = name.find('.')
        if search != -1:
            items = name.split('.')
            msg = "AutomaticDirichletConditionsProcess"
            msg += 'Please, specify the sub model part with the relative name\n'
            msg += '\tWrite : \"' + items[-1] + '\"\n'
            msg += '\tinstead of : \"' +  name + '\"\n'
            param.SetString(name)
            # The parameters are passed by value, we shall throw an error
            raise Exception(msg)
