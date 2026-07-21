import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MPMInletProcess(model, settings["Parameters"])

class MPMInletProcess(KratosMultiphysics.Process):
    def __init__(self, Model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        KratosMultiphysics.Process.__init__(self)

        # "modulus" must be a number
        if settings.Has("modulus") and not settings["modulus"].IsNumber():
            raise Exception('Parameter "modulus" must be a number')

        # "direction" must be a vector with three components
        if settings.Has("direction"):
            if not settings["direction"].IsVector() or settings["direction"].GetVector().Size() != 3:
                raise Exception('Parameter "direction" must be a vector of length 3')

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part_name     = settings["model_part_name"].GetString()
        self.model_part          = Model[self.model_part_name]
        self.assigned_velocity   = settings["velocity"].GetVector()
        self.gravity             = settings["gravity"].GetVector()
        self.height              = settings["height"].GetDouble()
        self.density             = settings["density"].GetDouble()
        self.vertical_distance   = settings["vertical_distance"].GetDouble()
        self.horizontal_distance = settings["horizontal_distance"].GetDouble()
        self.delta_time          = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        self.interval = KratosMultiphysics.IntervalUtility(settings)


    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name"       : "please_specify_model_part_name",
                "height"                : {},
                "density"               : 1.0,
                "vertical_distance"     : {},
                "horizontal_distance"   : {},
                "velocity"              : [0.0, 0.0, 0.0],
                "gravity"               : [0.0, 0.0, 0.0],
                "interval"              : [0.0, 1e30],
            }
            """)

    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):
            # Assign current velocities and displacement at the active nodes of the inlet
            self.step_is_active = True
            is_fixed = True
            self.variable_utils.ApplyFixity(KratosMultiphysics.VELOCITY    , is_fixed, self.model_part.Nodes)
            self.variable_utils.ApplyFixity(KratosMultiphysics.DISPLACEMENT, is_fixed, self.model_part.Nodes)

            self.variable_utils.SetVariable(KratosMultiphysics.VELOCITY    , self.assigned_velocity                  , self.model_part.Nodes)
            self.variable_utils.SetVariable(KratosMultiphysics.DISPLACEMENT, self.assigned_velocity * self.delta_time, self.model_part.Nodes)

            if (self._IsInletStep()):
                # Generate material points at the inlet

                # Assign (previous) velocities at the generated material points
                # Assign volume and mass to material points
                mp_volume = self.accumulated_volume / self.number_of_mp_to_be_generated
                mp_mass = mp_volume * self.density
                # for element in (vector of pointer of elements that was generated in this time step):
                    element.SetValuesOnIntegrationPoints(KratosMPM.MP_VELOCITY           ,[self.assigned_velocity],self.model_part.ProcessInfo)
                    element.SetValuesOnIntegrationPoints(KratosMPM.MP_VOLUME             ,[mp_volume]             ,self.model_part.ProcessInfo)
                    element.SetValuesOnIntegrationPoints(KratosMPM.MP_MASS               ,[mp_mass]               ,self.model_part.ProcessInfo)
                    element.SetValuesOnIntegrationPoints(KratosMPM.MP_VOLUME_ACCELERATION,[self.gravity]          ,self.model_part.ProcessInfo)
                #
            else:
                # if not inlet step, accumulate the volume and mass to be assigned to during the inlet step


        pass

    def ExecuteFinalizeSolutionStep(self):
        """This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.step_is_active:
            # Here we free all of the nodes in the mesh
            fixity_status  = False
            self.variable_utils.ApplyFixity(self.variable, fixity_status, self.model_part.Nodes)

        self.step_is_active = False

    def _IsInletStep(self):
        # Determine if this is the step of which material points will be generated at the inlet
        # It is determined by the spacing between material point given by the user
        # If the accumulated velocity * dt is bigger than the spacing, returns true
        # else, accumulate velocity * dt
        self.assigned_velocity *
        if ():
        pass