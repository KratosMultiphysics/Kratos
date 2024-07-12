import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def Factory(Settings, Model):
    if not isinstance(Settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string.")
    return ApplyWallLawProcess(Model, Settings["Parameters"])

class ApplyWallLawProcess(KratosMultiphysics.Process):

    @classmethod
    def _GetHelperClassSuffix(self):
        return "_monolithic_helper"

    class navier_slip_monolithic_helper():
        @classmethod
        def GetConditionName(cls):
            return "NavierStokesNavierSlipWallCondition"

        @classmethod
        def GetWallModelDefaultSettings(cls):
            wall_model_default_settings = KratosMultiphysics.Parameters("""{
                "slip_length" : 0.0
            }""")
            return wall_model_default_settings

        @classmethod
        def InitializeWallModelValues(cls, ModelPart, Settings):
            # Assign slip length to nodes
            slip_length = Settings["slip_length"].GetDouble()
            if not slip_length > 1.0e-12:
                # If provided value is zero, warn the user to set the SLIP_LENGTH manually
                warn_msg = "'slip_length' value found is zero. Positive non-zero value is expected to be set for 'SLIP_LENGTH' at the non-historical nodal database."
                KratosMultiphysics.Logger.PrintWarning(warn_msg)
            else:
                # If there is a provided value, assign it to the nodes
                for condition in ModelPart.Conditions:
                    for node in condition.GetNodes():
                        node.SetValue(KratosCFD.SLIP_LENGTH, slip_length)
                ModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(KratosCFD.SLIP_LENGTH)

    class linear_log_monolithic_helper():
        @classmethod
        def GetConditionName(self):
            return "NavierStokesLinearLogWallCondition"

        @classmethod
        def GetWallModelDefaultSettings(cls):
            wall_model_default_settings = KratosMultiphysics.Parameters("""{
                "y_wall" : 0.0
            }""")
            return wall_model_default_settings

        @classmethod
        def InitializeWallModelValues(cls, ModelPart, Settings):
            # Assign slip length to nodes
            y_wall = Settings["y_wall"].GetDouble()
            if not y_wall > 1.0e-12:
                # If provided value is zero, warn the user to set the Y_WALL manually
                warn_msg = "'y_wall' value found is zero. Positive non-zero value is expected to be set for 'Y_WALL' at the non-historical nodal database."
                KratosMultiphysics.Logger.PrintWarning(warn_msg)
            else:
                # If there is a provided value, assign it to the nodes
                for condition in ModelPart.Conditions:
                    for node in condition.GetNodes():
                        node.SetValue(KratosCFD.Y_WALL, y_wall)
                ModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(KratosCFD.Y_WALL)

    def __init__(self, Model, Settings):
        # Call base class constructor
        KratosMultiphysics.Process.__init__(self)

        # Validate input settings
        Settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        # Check input data
        if not Settings["model_part_name"].GetString():
            raise Exception("Empty 'model_part_name'. Provide a valid one.")
        if not Settings["wall_model_name"].GetString():
            err_msg = "Empty 'wall_model_name'. Available options are:"
            for wall_model_name in self.__GetAvailableModelsList():
                err_msg += f"\n\t-'{wall_model_name}'"
            raise Exception(err_msg)

        # Set the wall model helper class (see above)
        wall_model_name = Settings["wall_model_name"].GetString()
        self.wall_model_helper = self.__CreateWallModelHelperInstance(wall_model_name, self._GetHelperClassSuffix())

        # Validate the wall model settings with the corresponding wall model defaults
        Settings["wall_model_settings"].ValidateAndAssignDefaults(self.wall_model_helper.GetWallModelDefaultSettings())

        # Save auxiliary member variables
        self.model = Model
        self.settings = Settings

    def GetDefaultParameters(self):
        default_settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "",
            "wall_model_name" : "",
            "wall_model_settings" : {},
            "calculate_normals_at_each_step" : false
        }""")
        return default_settings

    def ExecuteInitialize(self):
        # Get data from settings
        model_part =  self.model[self.settings["model_part_name"].GetString()]

        # Get the condition registry name
        # Note that in here we are assuming a unique condition geometry in model part
        pts_num = 0
        for condition in model_part.Conditions:
            pts_num = condition.GetGeometry().PointsNumber()
            break
        pts_num = model_part.GetCommunicator().GetDataCommunicator().MaxAll(pts_num)
        domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        cond_reg_name = f"{self.wall_model_helper.GetConditionName()}{domain_size}D{pts_num}N"

        # Substitute current conditions by the corresponding ones implementing the wall model
        replace_settings = KratosMultiphysics.Parameters("""{
            "element_name" : "",
            "condition_name" : ""
        }""")
        replace_settings["condition_name"].SetString(cond_reg_name)
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(model_part, replace_settings).Execute()

        # Set the SLIP flag in the wall model part nodes and calculate the nodal normals
        # This is required in order to impose the no penetration condition by the rotating schemes
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(model_part, domain_size) #FIXME: This may interact with the nodal normals of slip boundaries (e.g. floor-wall in a 3D channel)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, model_part.Nodes)

        # Set the WALL flag in the wall model part new conditions
        # This is required in order to add the wall model contribution within the condition implementation
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.WALL, True, model_part.Conditions)

        # Initialize the corresponding wall model values using the wall model helper class (see above)
        self.wall_model_helper.InitializeWallModelValues(model_part, self.settings["wall_model_settings"])

    def ExecuteInitializeSolutionStep(self):
        # If required (e.g. moving boundaries) recalculate the nodal normals
        if self.settings["calculate_normals_at_each_step"].GetBool():
            model_part = self.model.GetSubModelPart(self.settings["model_part_name"].GetString())
            domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(model_part, domain_size) #FIXME: This may interact with the nodal normals of slip boundaries (e.g. floor-wall in a 3D channel)

    def __GetAvailableModelsList(self):
        available_models = []
        for attribute_name in dir(self):
            if self._GetHelperClassSuffix() in attribute_name:
                available_models.append(attribute_name.split(self._GetHelperClassSuffix())[0])
                # available_models.append(attribute_name.removesuffix(self._GetHelperClassSuffix())) #TODO: Use this after Python 3.9
        return available_models

    def __CreateWallModelHelperInstance(self, WallModelName, HelperSuffix):
        helper_class_name = f"{WallModelName}{HelperSuffix}"
        if not hasattr(self, helper_class_name):
            err_msg = f"Wrong wall model helper class name '{helper_class_name}'. Check your implementation."
            raise Exception(err_msg)
        return getattr(self, helper_class_name)()