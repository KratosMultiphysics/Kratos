import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignGravityToMaterialPointProcess(model, settings["Parameters"])

class AssignGravityToMaterialPointProcess(KratosMultiphysics.Process):
    def __init__(self, Model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        KratosMultiphysics.Process.__init__(self)

        # "modulus" must be a number
        if settings.Has("modulus") and not settings["modulus"].IsNumber():
            raise Exception('Parameter "modulus" must be a number')

        # "direction" must be a vector with three components
        if settings.Has("direction"):
            if not settings["direction"].IsVector() or settings["direction"].GetVector().Size() != 3:
                raise Exception('Parameter "direction" must be a vector of length 3')

        # "variable_name" is not requires: gravity is assigned
        # to variables "MP_ACCELERATION" and "MP_VOLUME_ACCELERATION"
        if settings.Has("variable_name"):
            settings.RemoveValue("variable_name")
            wrn_msg  = 'Parameter "variable_name" has been removed and will be ignored: '
            wrn_msg += 'gravity is assigned to "MP_VOLUME_ACCELERATION" and "MP_ACCELERATION".'
            KratosMultiphysics.Logger.PrintWarning("AssignGravityToMaterialPointProcess", wrn_msg)

        # This parameter was not used and is removed
        if settings.Has("constrained"):
            settings.RemoveValue("constrained")
            wrn_msg = 'Parameter "constrained" has been removed and will be ignored.'
            KratosMultiphysics.Logger.PrintWarning("AssignGravityToMaterialPointProcess", wrn_msg)

        # This parameter was not used and is removed
        if settings.Has("local_axes"):
            settings.RemoveValue("local_axes")
            wrn_msg = 'Parameter "local_axes" has been removed and will be ignored.'
            KratosMultiphysics.Logger.PrintWarning("AssignGravityToMaterialPointProcess", wrn_msg)

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.modulus = settings["modulus"].GetDouble()
        self.gravity_direction = settings["direction"].GetVector()
        self.gravity_acceleration = self.modulus * self.gravity_direction

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

        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            # Assign gravity to MP after solver.Initialize() - only apply once at the beginning!
            for element in self.model_part.Elements:
                element.SetValuesOnIntegrationPoints(KratosMPM.MP_VOLUME_ACCELERATION,[self.gravity_acceleration],self.model_part.ProcessInfo)
                element.SetValuesOnIntegrationPoints(KratosMPM.MP_ACCELERATION,[self.gravity_acceleration],self.model_part.ProcessInfo)
