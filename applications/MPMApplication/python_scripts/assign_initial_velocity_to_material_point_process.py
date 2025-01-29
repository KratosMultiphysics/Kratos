import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignInitialVelocityToMaterialPointProcess(model, settings["Parameters"])

class AssignInitialVelocityToMaterialPointProcess(KratosMultiphysics.Process):
    def __init__(self, Model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        KratosMultiphysics.Process.__init__(self)

        # "modulus" must be a number
        if settings.Has("modulus") and not settings["modulus"].IsNumber():
            raise Exception('Parameter "modulus" must be a number')

        # "direction" must be a vector with three components
        if settings.Has("direction"):
            if not settings["direction"].IsVector() or settings["direction"].GetVector().Size() != 3:
                raise Exception('Parameter "direction" must be a vector of length 3')

        # "variable_name" is not requires: velocity is assigned to variable "MP_VELOCITY"
        if settings.Has("variable_name"):
            settings.RemoveValue("variable_name")
            wrn_msg  = 'Parameter "variable_name" has been removed and will be ignored: '
            wrn_msg += 'initial velocity is assigned to "MP_VELOCITY".'
            KratosMultiphysics.Logger.PrintWarning("AssignInitialVelocityToMaterialPointProcess", wrn_msg)

        # This parameter was not used and is removed
        if settings.Has("constrained"):
            settings.RemoveValue("constrained")
            wrn_msg = 'Parameter "constrained" has been removed and will be ignored.'
            KratosMultiphysics.Logger.PrintWarning("AssignInitialVelocityToMaterialPointProcess", wrn_msg)

        # This parameter was not used and is removed
        if settings.Has("local_axes"):
            settings.RemoveValue("local_axes")
            wrn_msg = 'Parameter "local_axes" has been removed and will be ignored.'
            KratosMultiphysics.Logger.PrintWarning("AssignInitialVelocityToMaterialPointProcess", wrn_msg)

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        # Get updated model_part
        self.model = Model
        model_part_name = settings["model_part_name"].GetString()
        if (model_part_name.startswith('Initial_MPM_Material.')):
            model_part_name = model_part_name.replace('Initial_MPM_Material.','')
        self.mpm_material_model_part_name = "MPM_Material." + model_part_name
        # The actual initial velocity application occurs after the submodelpart is
        # transferred from the initial MPM material to the MPM material in the MaterialPointGeneratorUtility.
        # Therefore we change the prefix from initial MPM material
        # to MPM material.

        # Default settings
        self.modulus = settings["modulus"].GetDouble()
        self.velocity_direction = settings["direction"].GetVector()
        self.velocity = self.modulus * self.velocity_direction

    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name"      : "please_specify_model_part_name",
                "modulus"              : 1.0,
                "direction"            : [0.0, 0.0, 0.0]
            }
            """)

    def ExecuteBeforeSolutionLoop(self):
        # Assign velocity to MP after solver.Initialize() - only apply once at the beginning!
        model_part = self.model[self.mpm_material_model_part_name]
        # the model part is identified here, AFTER it has been transferred to the MPM_material part!
        if not model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            for element in model_part.Elements:
                element.SetValuesOnIntegrationPoints(KratosMPM.MP_VELOCITY,[self.velocity],model_part.ProcessInfo)
