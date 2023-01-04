import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def Factory(Settings, Model):
    if not isinstance(Settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string.")
    return ApplyWallLawProcess(Model, Settings["Parameters"])

class ApplyWallLawProcess(KratosMultiphysics.Process):

    __wall_model_condition_name_map = {
        "navier_slip" : "NavierStokesNavierSlipWallCondition",
        "linear_log" : "NavierStokesLinearLogWallCondition"
    }

    class navier_slip_helper():
        @classmethod
        def GetWallModelDefaultSettings(cls):
            wall_model_default_settings = KratosMultiphysics.Parameters("""{
                "slip_length" : 0.0
            }""")
            return wall_model_default_settings

        @classmethod
        def InitializeWallModelValues(cls, ModelPart, Settings):
            # Assign slip length to nodes
            slip_length = Settings["slip_lenght"].GetDouble()
            if not slip_length > 1.0e-12:
                # If provided value is zero, warn the user to set the SLIP_LENGTH manually
                warn_msg = f"'slip_length' value found is zero. Positive non-zero value is expected to be set for 'SLIP_LENGTH' at the non-historical nodal database."
                KratosMultiphysics.Logger.PrintWarning(warn_msg)
            else:
                # If there is a provided value, assign it to the nodes
                for condition in ModelPart.Conditions:
                    for node in condition.GetNodes:
                        node.SetValue(KratosCFD.SLIP_LENGTH, slip_length)
                ModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(KratosCFD.SLIP_LENGTH)

    class linear_log_helper():
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
                warn_msg = f"'y_wall' value found is zero. Positive non-zero value is expected to be set for 'Y_WALL' at the non-historical nodal database."
                KratosMultiphysics.Logger.PrintWarning(warn_msg)
            else:
                # If there is a provided value, assign it to the nodes
                for condition in ModelPart.Conditions:
                    for node in condition.GetNodes:
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
            for wall_model_name in self.__wall_model_condition_name_map.keys():
                err_msg += f"\n\t-'{wall_model_name}'"
            raise Exception(err_msg)

        # Set the wall model helper class (see above)
        wall_model_name = Settings["wall_model_name"].GetString()
        helper_class_name = f"{wall_model_name}_helper"
        if not hasattr(self, helper_class_name):
            err_msg = f"Wrong wall model helper class name '{helper_class_name}'. Expected in this case is '__{wall_model_name}_helper'. Check your implementation."
            raise Exception(err_msg)
        self.wall_model_helper = getattr(self, helper_class_name)()

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
        wall_model_name = self.settings["wall_model_name"].GetString()

        # Get the condition registry name
        # Note that in here we are assuming a unique condition geometry in model part
        pts_num = 0
        for condition in model_part.Conditions:
            pts_num = condition.PointsNumber()
            break
        pts_num = model_part.GetCommunicator().GetDataCommunicator().MaxAll(pts_num)
        domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        cond_aux_name = self.__wall_model_condition_name_map[wall_model_name]
        cond_reg_name = f"{cond_aux_name}{domain_size}D{pts_num}N"

        # Substitute current conditions by the corresponding ones implementing the wall model
        # 1. Save current condition data in a map such that Id : (geometry, properties)
        aux_cond_data_map = {}
        for condition in model_part.Conditions:
            condition.Set(KratosMultiphysics.TO_ERASE,True)
            aux_cond_data_map[condition.Id] = (condition.GetGeometry(), condition.Properties)
        # 2. Remove current conditions
        # Note that in here we need to assume that no submodelpart is hanging from this one as the hierarchical structure would be destroyed
        if model_part.NumberOfSubModelParts() != 0:
            raise Exception(f"Wall model part '{model_part.FullName}' has submodelparts.")
        KratosMultiphysics.AuxiliarModelPartUtilities(model_part).RemoveConditionsAndBelongings(KratosMultiphysics.TO_ERASE)
        # 3. Create new wall law conditions from the data map info
        for key, value in aux_cond_data_map.items():
            model_part.CreateNewCondition(cond_reg_name, key, value[0], value[1])

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