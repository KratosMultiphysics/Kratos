# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsBiomedicalApplication as KratosBio

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyParabolicInletProcess(model, settings["Parameters"])

class ApplyParabolicInletProcess(KratosMultiphysics.Process):
    """
    Class to compute parabolic inlet profile
    """

    def __init__(self, model, settings):
        # Base class constructor call
        KratosMultiphysics.Process.__init__(self)

        # Set default settings
        # Note the trick to allow both sacalar and function values of the parabola maximum value
        default_settings = self.GetDefaultParameters()
        if settings.Has("parabola_vertex_value"):
            if settings["parabola_vertex_value"].IsString():
                settings["parabola_vertex_value"].SetString("0.0")
        else:
            raise Exception("'max_value' not found. It needs to be user-provided.")

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # Check default settings
        settings.ValidateAndAssignDefaults(default_settings)

        # Check user-provided data
        if not settings["inlet_model_part_name"].GetString():
            raise ValueError("'inlet_model_part' not provided.")

        if not settings["wall_model_part_name"].GetString():
            raise ValueError("'wall_model_part' needs to be provided for 'parabolic' inlet distribution.")

        # Set the maximum value input
        self.max_value_is_numeric = False
        if settings["parabola_vertex_value"].IsNumber():
            self.max_value_is_numeric = True
            self.max_value = settings["parabola_vertex_value"].GetDouble()
        else:
            self.function_string = settings["parabola_vertex_value"].GetString()
            self.max_value_function = KratosMultiphysics.GenericFunctionUtility(self.function_string, settings["local_axes"])

        # Save model and settings containers
        self.model = model
        self.settings = settings

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters("""{
            "wall_model_part_name": "",
            "inlet_model_part_name": "",
            "parabola_vertex_value" : 0.0,
            "interval" : [0.0,"End"],
            "local_axes" : {}
        }""")
        return default_settings

    def ExecuteInitialize(self):
        # Get and check domain size
        inlet_model_part = self.model.GetModelPart(self.settings["inlet_model_part_name"].GetString())
        domain_size = inlet_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2,3]:
            raise ValueError(f"Wrong 'DOMAIN_SIZE' value {domain_size} in ProcessInfo container.")

        # Set the INLET flag in the inlet model part nodes and conditions
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INLET, True, inlet_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INLET, True, inlet_model_part.Conditions)

        # Compute the normal on the nodes of interest
        # Note that a custom normal variable is use to avoid interfering with the wall one
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplexNonHistorical(
            inlet_model_part,
            domain_size,
            KratosBio.INLET_NORMAL)

        # Calculate the distance to the wall
        # This will be used to calculate the inlet parabolic profile
        #FIXME: HERE WE SHOULD USE THE PARALLEL DISTANCE CALCULATOR
        wall_model_part = self.model.GetModelPart(self.settings["wall_model_part_name"].GetString())
        root_model_part = wall_model_part.GetRootModelPart()
        distance_settings = KratosMultiphysics.Parameters("""{
            "distance_variable" : "WALL_DISTANCE",
			"distance_database" : "nodal_non_historical"
        }""")
        distance_process = self._ReturnDistanceProcessPrototype(domain_size)(
            root_model_part,
            wall_model_part,
            distance_settings
        )
        distance_process.Execute()

    def ExecuteBeforeSolutionLoop(self):
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        inlet_model_part = self.model.GetModelPart(self.settings["inlet_model_part_name"].GetString())
        current_time = inlet_model_part.ProcessInfo[KratosMultiphysics.TIME]
        if self.interval.IsInInterval(current_time):
            self.step_is_active = True
            if self.max_value_is_numeric:
                KratosBio.ParabolicProfileUtilities.ImposeParabolicInlet(inlet_model_part, self.max_value)
            else:
                KratosBio.ParabolicProfileUtilities.ImposeParabolicInlet(inlet_model_part, self.max_value_function)

    def ExecuteFinalizeSolutionStep(self):
        # Here we free all of the nodes in the inlet
        if self.step_is_active:
            inlet_model_part = self.model.GetModelPart(self.settings["inlet_model_part_name"].GetString())
            KratosBio.ParabolicProfileUtilities.FreeParabolicInlet(inlet_model_part)
        self.step_is_active = False

    @classmethod
    def _ReturnDistanceProcessPrototype(cls, domain_size):
        if domain_size == 2:
            return KratosMultiphysics.CalculateDistanceToSkinProcess2D
        else:
            return KratosMultiphysics.CalculateDistanceToSkinProcess3D
