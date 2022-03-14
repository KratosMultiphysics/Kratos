import KratosMultiphysics


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return DeceleratingWallProcess(model, settings["Parameters"])


class DeceleratingWallProcess(KratosMultiphysics.Process):
    """
    This process helps with convergence by slowly transforming an outlet into a
    wall. During the start of the simulation, the normal component of the
    momentum of the selected wall will have exponential decay. After its speed
    is 1000 times smaller than originally, it becomes a wall.

    The effect on the normal component of the momentum `f`, if only manipulated by
    this process is:
    ```
           { f0 * 1000^(-t/period)  if t < period
    f(t) = {
           { 0                      if t >= period
    ```
    If the momentum is being affected by any other variable, then the long-term
    effect is less predictable but each step the momentum is diminished by a
    factor of 1000^(-dt/period).

    Key parameters
    - Period: For a variable with otherwise constant momentum, it's the time it
    takes for the momentum to be 0.001 of the original. After this period, the
    model part is treated like a solid wall.
    """

    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        assert settings.Has("period")

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        modelpart_name = settings["model_part_name"].GetString()
        self.model_part = model[modelpart_name]

        self.period = settings["period"].GetDouble()
        self.recompute_normals = settings["recompute_normals"].GetBool()

        if self.period <= 0:
            raise ValueError("Parameter period must be strictly positive.")

        self.period_ended = False

    def GetDefaultParameters(self):
        return KratosMultiphysics.Parameters(
            """
            {
                "model_part_name" : "SPECIFY_MODEL_PART_NAME",
                "period"        : 0,
                "recompute_normals" : false
            }
            """
        )

    def ExecuteInitialize(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateNormals(self.model_part, True)

    def ExecuteFinalizeSolutionStep(self):
        if self.period_ended:
            return

        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if time > self.period:
            KratosMultiphysics.VariableUtils().SetFlag(
                KratosMultiphysics.SLIP, True, self.model_part.Nodes
            )
            self.period_ended = True
            return

        if self.recompute_normals:
            KratosMultiphysics.NormalCalculationUtils().CalculateNormals(self.model_part, True)

        dt = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        decay = (1000)**(-dt/self.period)

        for node in self.model_part.Nodes:
            mom = node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM)

            unit_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            unit_normal /= unit_normal.norm_2()

            mom_n = sum([m*n for m, n in zip(mom, unit_normal)])
            new_mom_n = mom_n * decay

            delta_mom_n = new_mom_n - mom_n

            new_mom = mom + delta_mom_n * unit_normal

            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM, new_mom)
