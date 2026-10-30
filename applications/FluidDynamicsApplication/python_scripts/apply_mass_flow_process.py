import numpy

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMassFlowProcess(Model, settings["Parameters"])


class ApplyMassFlowProcess(KratosMultiphysics.Process):
    """This process imposes a mass flow on a list of velocity boundaries.

    For each boundary the nodal VELOCITY is scaled such that the flow rate through its conditions
    (as computed by FluidAuxiliaryUtilities.CalculateFlowRate) matches the prescribed mass flow.
    As this is the flow rate appearing in the discrete mass conservation equation, prescribing
    mass flows that add up to zero yields a solvable problem when all the boundaries have the
    velocity imposed. Note that the match is exact for linear simplicial conditions.

    The "mass_flow" can be a number or a function of time, being positive when entering the domain.
    In 2D it is the mass flow per unit depth.

    Two velocity profiles are available for each boundary:
    - "uniform": the process sets (and fixes if "constrained") a uniform velocity along the average
      boundary normal. The velocity components already fixed (e.g. by a wall or a symmetry condition) are
      left untouched and their flow rate is subtracted from the prescribed one.
    - "rescale_existing": the process scales the velocity set by the previous processes.

    Hence, this process must be placed after the wall processes and after the processes setting
    the profile to be rescaled. Boundaries can share nodes (e.g. adjacent outlets) if these use the
    "uniform" profile and are "constrained": the shared nodes take the velocity of the first boundary
    in the list and the next ones treat them as fixed.
    """

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"   : "",
            "boundaries"        : [],
            "density"           : 0.0,
            "constrained"       : true,
            "check_balance"     : true,
            "balance_tolerance" : 1.0e-10,
            "interval"          : [0.0,"End"]
        }
        """)

        # If not provided, the density is taken from the properties of the fluid model part elements
        user_defined_density = settings.Has("density")

        settings.ValidateAndAssignDefaults(default_settings)

        if settings["model_part_name"].GetString() == "":
            raise Exception("Empty fluid model part name string. Set a valid model part name.")
        if settings["boundaries"].size() == 0:
            raise Exception("No boundaries are specified in 'boundaries'.")

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.constrained = settings["constrained"].GetBool()
        self.check_balance = settings["check_balance"].GetBool()
        self.balance_tolerance = settings["balance_tolerance"].GetDouble()
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        self.density = settings["density"].GetDouble() if user_defined_density else self.__GetDensityFromProperties()
        if self.density <= 0.0:
            raise Exception(f"DENSITY is {self.density}, positive number expected.")

        self.boundaries = [self.__CreateBoundary(Model, boundary_settings) for boundary_settings in settings["boundaries"].values()]

    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        if not self.interval.IsInInterval(current_time):
            return

        # Evaluate the mass flows and check that these add up to zero
        mass_flows = [boundary["mass_flow"](current_time) for boundary in self.boundaries]
        if self.check_balance:
            net_mass_flow = sum(mass_flows)
            max_mass_flow = max(abs(mass_flow) for mass_flow in mass_flows)
            if abs(net_mass_flow) > self.balance_tolerance * max_mass_flow:
                err_msg = f"Prescribed mass flows do not add up to zero at time {current_time} (net mass flow is {net_mass_flow}).\n"
                for boundary, mass_flow in zip(self.boundaries, mass_flows):
                    err_msg += f"\t- '{boundary['model_part'].FullName()}': {mass_flow}\n"
                err_msg += "Set 'check_balance' to false if this is intended (e.g. there is a pressure outlet)."
                raise Exception(err_msg)

        for boundary, mass_flow in zip(self.boundaries, mass_flows):
            self.__ImposeMassFlow(boundary, mass_flow)

    def ExecuteFinalizeSolutionStep(self):
        # Free the DOFs fixed by this process
        for boundary in self.boundaries:
            own_dofs = boundary["fixed_dofs"]
            if own_dofs is not None:
                fixity_ta = boundary["fixity_ta"]
                fixity_ta.CollectData()
                fixity = fixity_ta.data.copy()
                fixity[own_dofs] = False
                fixity_ta.data = fixity
                fixity_ta.StoreData()
                boundary["fixed_dofs"] = None

    def __CreateBoundary(self, Model, boundary_settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "",
            "mass_flow"       : 0.0,
            "profile"         : "uniform"
        }
        """)

        # Trick: allow "mass_flow" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if boundary_settings.Has("mass_flow"):
            if boundary_settings["mass_flow"].IsString():
                default_settings["mass_flow"].SetString("0.0")
        else:
            raise Exception("Please specify the 'mass_flow' of each boundary.")

        boundary_settings.ValidateAndAssignDefaults(default_settings)

        if boundary_settings["model_part_name"].GetString() == "":
            raise Exception("Empty boundary model part name string. Set a valid model part name.")
        boundary_model_part = Model[boundary_settings["model_part_name"].GetString()]

        profile = boundary_settings["profile"].GetString()
        if profile not in ["uniform", "rescale_existing"]:
            raise Exception(f"Wrong 'profile' in '{boundary_model_part.FullName()}'. Provided value is '{profile}'. Available options are 'uniform' and 'rescale_existing'.")

        if boundary_settings["mass_flow"].IsNumber():
            mass_flow_value = boundary_settings["mass_flow"].GetDouble()
            mass_flow = lambda time : mass_flow_value
        else:
            mass_flow_function = KratosMultiphysics.GenericFunctionUtility(boundary_settings["mass_flow"].GetString())
            if mass_flow_function.DependsOnSpace():
                raise Exception(f"The 'mass_flow' of '{boundary_model_part.FullName()}' can only depend on time.")
            mass_flow = lambda time : mass_flow_function.CallFunction(0.0, 0.0, 0.0, time, 0.0, 0.0, 0.0)

        boundary = {
            "model_part" : boundary_model_part,
            "mass_flow" : mass_flow,
            "profile" : profile,
            "area" : self.__CalculateArea(boundary_model_part),
            "velocity_ta" : KratosMultiphysics.TensorAdaptors.HistoricalVariableTensorAdaptor(boundary_model_part.Nodes, KratosMultiphysics.VELOCITY),
            "fixity_ta" : KratosMultiphysics.TensorAdaptors.FixityTensorAdaptor(boundary_model_part.Nodes, [KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.VELOCITY_Z]),
            "fixed_dofs" : None
        }

        if profile == "uniform":
            boundary["inwards_normal"] = self.__CalculateInwardsNormal(boundary_model_part)

        return boundary

    def __ImposeMassFlow(self, boundary, mass_flow):
        boundary_model_part = boundary["model_part"]
        velocity_ta = boundary["velocity_ta"]
        velocity_ta.CollectData()
        velocity = velocity_ta.data.copy()

        # Set the velocity profile to be scaled
        if boundary["profile"] == "uniform":
            # Note that the velocity components already fixed (e.g. those of a wall or a symmetry plane) are not modified
            # Their flow rate is computed with the free components set to zero and subtracted from the prescribed one
            fixity_ta = boundary["fixity_ta"]
            fixity_ta.CollectData()
            fixity = fixity_ta.data.copy()
            own_dofs = ~fixity
            velocity[own_dofs] = 0.0
            velocity_ta.data = velocity
            velocity_ta.StoreData()
            fixed_flow_rate = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRate(boundary_model_part)
            velocity = numpy.where(own_dofs, boundary["inwards_normal"], velocity)
            velocity_ta.data = velocity
            velocity_ta.StoreData()
            if self.constrained:
                fixity[own_dofs] = True
                fixity_ta.data = fixity
                fixity_ta.StoreData()
                boundary["fixed_dofs"] = own_dofs
        else:
            own_dofs = slice(None)
            fixed_flow_rate = 0.0

        # Scale the profile to match the prescribed mass flow
        # Note that the flow rate is computed with the outwards normal (i.e. it is negative if entering the domain)
        if mass_flow == 0.0 and fixed_flow_rate == 0.0:
            scale_factor = 0.0
        else:
            profile_flow_rate = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRate(boundary_model_part) - fixed_flow_rate
            # Note that the reference flow rate to check if the profile one is zero is the max. velocity times the boundary area
            max_velocity = float(abs(velocity).max()) if velocity.size != 0 else 0.0
            max_velocity = boundary_model_part.GetCommunicator().GetDataCommunicator().MaxAll(max_velocity)
            if abs(profile_flow_rate) <= 1.0e-12 * max_velocity * boundary["area"]:
                err_msg = f"The '{boundary['profile']}' profile in '{boundary_model_part.FullName()}' has zero flow rate, so it cannot be scaled to the prescribed mass flow.\n"
                err_msg += "Check that the velocity profile is set (or that not all the nodes are fixed) before this process is executed."
                raise Exception(err_msg)
            scale_factor = (-mass_flow / self.density - fixed_flow_rate) / profile_flow_rate

        velocity[own_dofs] *= scale_factor
        velocity_ta.data = velocity
        velocity_ta.StoreData()
        boundary_model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.VELOCITY)

    def __GetDensityFromProperties(self):
        density = 0.0
        for element in self.model_part.GetCommunicator().LocalMesh().Elements:
            density = element.Properties.GetValue(KratosMultiphysics.DENSITY)
            break
        return self.model_part.GetCommunicator().GetDataCommunicator().MaxAll(density)

    @staticmethod
    def __CalculateInwardsNormal(boundary_model_part):
        # Compute the condition and nodes area normals and save the values in the NORMAL variable
        # Also note that in the nodes case, the normal values are stored in the historical database
        enforce_generic_algorithm = False
        KratosMultiphysics.NormalCalculationUtils().CalculateNormals(
            boundary_model_part,
            enforce_generic_algorithm,
            KratosMultiphysics.NORMAL)

        # Compute the average conditions normal in the boundary
        area_normal = KratosMultiphysics.VariableUtils().SumConditionVectorVariable(KratosMultiphysics.NORMAL, boundary_model_part)
        area = area_normal.norm_2()
        if area < 1.0e-12:
            err_msg = f"Area normal norm is close to 0 in '{boundary_model_part.FullName()}'. This might be caused because:\n"
            err_msg += "\t- the provided model part has no conditions\n"
            err_msg += "\t- the provided model part is a closed surface"
            raise Exception(err_msg)

        # Note that it is assumed that the NormalCalculationUtils().CalculateNormals returns the outwards normal vector
        return [-area_normal[i] / area for i in range(3)]

    @staticmethod
    def __CalculateArea(boundary_model_part):
        area = 0.0
        for condition in boundary_model_part.GetCommunicator().LocalMesh().Conditions:
            area += condition.GetGeometry().DomainSize()
        return boundary_model_part.GetCommunicator().GetDataCommunicator().SumAll(area)
