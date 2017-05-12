import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyWallFunctionProcess(Model, settings["Parameters"])

import apply_slip_process

class ApplyWallFunctionProcess(apply_slip_process.ApplySlipProcess):
    def __init__(self, Model, settings):

        default_parameters = KratosMultiphysics.Parameters("""
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "avoid_recomputing_normals": false,
                "wall_model": "fs_werner_wengle"
            }""")

        settings.ValidateAndAssignDefaults(default_parameters)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.wall_model = settings["wall_model"].GetString()

        # initialize Y_WALL=0.0 and IS_STRUCTURE=1.0 on nodes
        slip_settings = KratosMultiphysics.Parameters(settings)
        slip_settings.RemoveValue("wall_model")
        super(ApplyWallFunctionProcess, self).__init__(Model, slip_settings)
        # get wall model settings
        domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if self.wall_model == "fs_werner_wengle":
            condition_name = "FSWernerWengleWallCondition" + str(domain_size) + "D"
        elif self.wall_model == "fs_generalized":
            condition_name = "FSGeneralizedWallCondition" + str(domain_size) + "D"
        else:
            raise Exception("invalid wall_model: " + self.wall_model)

        self.replacement_settings = KratosMultiphysics.Parameters("{}")
        self.replacement_settings.AddEmptyValue("condition_name")
        self.replacement_settings["condition_name"].SetString(condition_name)

    def ExecuteInitialize(self):
        # replace conditions of this (sub) model part with wall function
        KratosMultiphysics.ReplaceConditionsProcess(self.model_part, self.replacement_settings).Execute()
        # set flags to activate wall function
        for cond in self.model_part.Conditions:
            cond.SetValue(KratosMultiphysics.IS_STRUCTURE,1.0)
            cond.SetValue(KratosMultiphysics.Y_WALL, 1.0) # Y_WALL > 0 to activate wall model
        for node in self.model_part.Nodes:
            node.SetValue(KratosMultiphysics.Y_WALL, 1.0) # Y_WALL > 0 to activate wall model
        # used by conditions to find parent element
        KratosMultiphysics.FindNodalNeighboursProcess(self.model_part.GetRootModelPart(),10,10).Execute()
        # needed by some wall functions
        normal_util = KratosMultiphysics.NormalCalculationUtils()
        normal_util.CalculateOnSimplex(self.model_part.GetRootModelPart(),self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
