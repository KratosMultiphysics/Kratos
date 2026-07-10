import KratosMultiphysics


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CompressibleSlipWallProcess(model, settings["Parameters"])


class CompressibleSlipWallProcess(KratosMultiphysics.Process):
    """
    This process helps with convergence by slowly transforming an outlet into a
    wall. This process has four phases:

    Stages
    ------
    - `AWAITING`: Nothing is done before the interval starts.
    - `RAMP_UP`: The time after interval start. Duration is rampup_time. The
    wall lets a proportion of the fluid through.
    - `STEADY` state: Solid free-slip wall. Lasts until the end of the interval.
    - `FINISHED`: Nothing is done after the interval ends.

    RAMP-UP
    -------
    During this stage, the normal component of the momentum or velocity of the
    selected wall will have exponential decay. After its speed is A times
    smaller than originally, it exits ramp-up stage.

    This parameter A is hardcoded as `decay_constant`.

    The effect on the normal component of the momentum/velocity `f`, if only
    manipulated by this process is:
    ```Python
           { f0 * A^(-t/period)  if t < period (RAMP-UP)
    f(t) = {
           { 0                   if t >= period (STEADY STATE)
    ```
    If the momentum/velocity is being affected by any other variable, then the
    long-term effect is less predictable but each step the momentum/velocity is
    diminished by a factor of A^(-dt/period).

    STEADY STATE
    ------------
    The node is given the SLIP flag, which signals the strategy that it is in
    charge of applying the boundary condition.

    Key parameters
    --------------
    - rampup_time: Duration of the ramp-up stage.
    - interval:  The time in this process runs (RAMP_UP + STEADY).
    - recompute_normals: Enable computing the normals every step.
    - variable: The variable to control, usually MOMENTUM of VELOCITY
    """

    decay_constant = 1000

    # Stages
    AWAITING = 0
    RAMP_UP = 1
    STEADY = 2
    FINISHED = 3

    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        self.interval = KratosMultiphysics.IntervalUtility(settings)

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        modelpart_name = settings["model_part_name"].GetString()

        self.model_part = model[modelpart_name]
        self.rampup_time = settings["rampup_time"].GetDouble()
        self.rampup_enabled = True
        self.recompute_normals = settings["recompute_normals"].GetBool()
        self._stage = self.AWAITING
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())

    def GetDefaultParameters(self):
        return KratosMultiphysics.Parameters(
            """
            {
                "model_part_name"   : "SPECIFY_MODEL_PART_NAME",
                "interval"   : [0, "End"],
                "rampup_time"   : 0,
                "recompute_normals" : false,
                "variable_name" : "MOMENTUM"
            }
            """
        )

    def Check(self):
        if self.rampup_time <= 0:
            raise ValueError("Parameter rampup_time must be positive or zero.")

        full_duration = self.interval.GetIntervalEnd() - self.interval.GetIntervalBegin()
        if self.rampup_time > full_duration:
            msg = "Parameter rampup_time must smaller than the interval ({0} < {1})."\
                .format(self.rampup_time, full_duration)
            raise ValueError(msg)

        if not self.model_part.HasNodalSolutionStepVariable(self.variable):
            msg = "The variable must be in the historical database. Provided variable {} is not.".format(self.variable.Name())
            raise RuntimeError(msg)

        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(self.variable.Name())
        normal_type = KratosMultiphysics.KratosGlobals.GetVariableType(KratosMultiphysics.NORMAL.Name())
        if variable_type != normal_type:
            msg = "Variable must be a vector with the same size as NORMAL. Provided variable {} is not."\
                .format(self.variable.Name())
            raise TypeError(msg)

    def ExecuteInitialize(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateNormals(self.model_part, True)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, False, self.model_part.Nodes)

        if not self.rampup_time > 0.0:
            self.rampup_enabled = False

    def ExecuteInitializeSolutionStep(self):
        self._UpdateStage()

        if self.recompute_normals and self.stage in [self.RAMP_UP, self.STEADY]:
            KratosMultiphysics.NormalCalculationUtils().CalculateNormals(self.model_part, True)

    def ExecuteFinalizeSolutionStep(self):
        if self._stage != self.RAMP_UP:
            return

        dt = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        decay = (self.decay_constant)**(-dt/self.rampup_time)
        div_by_zero_tolerance = 1e-8

        for node in self.model_part.Nodes:
            unit_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            norm2 = unit_normal.norm_2()
            if norm2 < div_by_zero_tolerance:
                continue
            unit_normal /= norm2

            vec = node.GetSolutionStepValue(self.variable)
            vec_n = sum([m*n for m, n in zip(vec, unit_normal)])
            new_vec_n = vec_n * decay

            delta_vec_n = new_vec_n - vec_n
            new_vec = vec + delta_vec_n * unit_normal

            node.SetSolutionStepValue(self.variable, new_vec)

    def _GetStage(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        if time <= self.interval.GetIntervalBegin():
            return self.AWAITING

        rampup_stop = self.interval.GetIntervalBegin() + self.rampup_time
        if self.rampup_enabled and time <= rampup_stop:
            return self.RAMP_UP

        elif time <= self.interval.GetIntervalEnd():
            return self.STEADY
        return self.FINISHED

    def _UpdateStage(self):
        new_stage = self._GetStage()

        if new_stage == self._stage:
            return

        if new_stage == self.STEADY:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, self.model_part.Nodes)
        elif new_stage == self.FINISHED:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, False, self.model_part.Nodes)

        self._stage = new_stage
