import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

from KratosMultiphysics import assign_vector_by_direction_process

def Factory(Settings, Model):
    if not isinstance(Settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string.")
    return ApplyWallLawProcess(Model, Settings["Parameters"])

class ApplyWallLawProcess(KratosMultiphysics.Process):

    __wall_model_condition_name_map = {
        "navier_slip" : "NavierStokesNavierSlipWallCondition",
        "linear_log" : "NavierStokesLinearLogWallCondition"
    }

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

        # Save member variables
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
        # Get wall model part
        model_part =  self.Model[self.settings["model_part_name"].GetString()]

        # Get the condition registry name
        # Note that in here we are assuming a unique condition geometry in model part
        pts_num = 0
        for condition in model_part.Conditions:
            pts_num = condition.PointsNumber()
            break
        model_part.GetCommunicator().GetDataCommunicator().MaxAll(pts_num)
        domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        cond_aux_name = self.__wall_model_condition_name_map[self.settings["wall_model_name"].GetString()]
        cond_reg_name = f"{cond_aux_name}{domain_size}D{pts_num}N"

        # Substitute current conditions by the corresponding ones implementing the wall model
        # 1. Save current condition data in a map such that Id : (geometry, properties)
        aux_cond_data_map = {}
        for condition in model_part.Conditions:
            condition.Set(KratosMultiphysics.TO_ERASE,True)
            aux_cond_data_map[condition.Id] = (condition.GetGeometry(), condition.Properties)
        # 2. Remove current conditions
        # Note that in here we need to assume that no submodelpart is hanging from this one as the hierarchical structure would be destroyed
        if model_part.NumberOfSubModelParts != 0:
            raise Exception(f"Wall model part '{model_part.FullName}' has submodelparts.")
        KratosMultiphysics.AuxiliarModelPartUtilities(model_part).RemoveConditionsAndBelongings(KratosMultiphysics.TO_ERASE)
        # 3. Create new wall law conditions from the data map info
        for key, value in aux_cond_data_map:
            model_part.CreateNewCondition(cond_reg_name, key, value[0], value[1])

        # Set the SLIP flag in the wall model part nodes and calculate the nodal normals
        # This is required in order to impose the no penetration condition by the rotating schemes
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(model_part, domain_size) #FIXME: This may interact with the nodal normals of slip boundaries (e.g. floor-wall in a 3D channel)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, model_part.Nodes)

        # Set the WALL flag in the wall model part new conditions
        # This is required in order to add the wall model contribution within the condition implementation
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.WALL, True, model_part.Conditions)

        # Set the wall law required data
        #TODO: Create helper classes similar to what we do in the elements

    def ExecuteInitializeSolutionStep(self):
        # If required (e.g. moving boundaries) recalculate the nodal normals
        if self.settings["calculate_normals_at_each_step"].GetBool():
            self.ExecuteInitialize()
