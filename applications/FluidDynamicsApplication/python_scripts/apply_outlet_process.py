import math
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyOutletProcess(Model, settings["Parameters"])


class ApplyOutletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "mesh_id"            : 0,
            "model_part_name"    : "",
            "variable_name"      : "PRESSURE",
            "constrained"        : true,
            "value"              : 0.0,
            "interval"           : [0.0,"End"],
            "hydrostatic_outlet" : false,
            "h_top"              : 0.0
        }
        """)

        # Trick: allows "value" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if(settings.Has("value")):
            if(settings["value"].IsString()):
                default_settings["value"].SetString("0.0")

        settings.ValidateAndAssignDefaults(default_settings)

        # Set a Kratos parameters suitable for the core processes to set the PRESSURE
        pres_settings = settings.Clone()
        pres_settings.RemoveValue("hydrostatic_outlet")
        pres_settings.RemoveValue("h_top")

        # Create a copy of the PRESSURE settings to set the EXTERNAL_PRESSURE
        ext_pres_settings = pres_settings.Clone()
        ext_pres_settings["constrained"].SetBool(False)
        ext_pres_settings["variable_name"].SetString("EXTERNAL_PRESSURE")

        # Check the core processes input data
        if (pres_settings["model_part_name"].GetString() == ""):
            raise Exception("Empty outlet pressure model part name. Set a valid model part name.")
        elif (ext_pres_settings["model_part_name"].GetString() == ""):
            raise Exception("Empty outlet external pressure model part name. Set a valid model part name.")
        elif (pres_settings["variable_name"].GetString() != "PRESSURE"):
            raise Exception("Outlet pressure settings variable_name is not PRESSURE.")
        elif (ext_pres_settings["variable_name"].GetString() != "EXTERNAL_PRESSURE"):
            raise Exception("Outlet external pressure settings variable_name is not EXTERNAL_PRESSURE.")
        elif (pres_settings["value"].IsString()):
            if (pres_settings["value"].GetString == ""):
                raise Exception("Outlet pressure function sting is empty.")
        elif (ext_pres_settings["value"].IsString()):
            if (ext_pres_settings["value"].GetString == ""):
                raise Exception("Outlet external pressure function sting is empty.")

        self.hydrostatic_outlet = settings["hydrostatic_outlet"].GetBool()
        self.h_top = settings["h_top"].GetDouble()

                # Set the OUTLET flag in the outlet model part nodes and conditions
        self.outlet_model_part = Model[pres_settings["model_part_name"].GetString()]
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.outlet_model_part, self.outlet_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
        for node in self.outlet_model_part.Nodes:
            node.Set(KratosMultiphysics.OUTLET, True)
        for condition in self.outlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.OUTLET, True)

        # Construct the base process AssignValueProcess
        import assign_scalar_variable_process
        self.aux_pressure_process = assign_scalar_variable_process.AssignScalarVariableProcess(Model, pres_settings)
        self.aux_external_pressure_process = assign_scalar_variable_process.AssignScalarVariableProcess(Model, ext_pres_settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        self.aux_pressure_process.ExecuteInitializeSolutionStep()
        self.aux_external_pressure_process.ExecuteInitializeSolutionStep()

        # If considered, add the hydrostatic component to the outlet pressure
        if (self.hydrostatic_outlet):
            self._AddOutletHydrostaticComponent()

        # Compute the outlet average velocity
        self._ComputeOutletCharacteristicVelocity()


    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        self.aux_pressure_process.ExecuteFinalizeSolutionStep()
        self.aux_external_pressure_process.ExecuteFinalizeSolutionStep()


    # Private methods section
    def _AddOutletHydrostaticComponent(self):
        # Initialize body force value (avoid segfault in MPI if the local mesh has no outlet nodes)
        body_force = KratosMultiphysics.Vector(3)
        body_force[0] = 0.0
        body_force[1] = 0.0
        body_force[2] = 0.0

        # Get the body force value
        for node in self.outlet_model_part.Nodes:
            body_force = node.GetSolutionStepValue(KratosMultiphysics.BODY_FORCE, 0)
            break

        # Compute the body force unit normal vector
        body_force_norm = math.sqrt(body_force[0]*body_force[0] + body_force[1]*body_force[1] + body_force[2]*body_force[2])    # Body force norm
        body_force_dir = (1/(body_force_norm+1e-10))*body_force                                                                 # Body force unit director vector

        # Compute the minimum body force projection value (reference value)
        min_proj = 0.0
        for node in self.outlet_model_part.Nodes:
            body_force_proj = body_force_dir[0]*node.X + body_force_dir[1]*node.Y + body_force_dir[2]*node.Z    # Iteration node body force projection
            min_proj = min(min_proj, body_force_proj)

        # Add the hydrostatic component to the current PRESSURE and EXTERNAL_PRESSURE values
        for node in self.outlet_model_part.Nodes:
            body_force_proj = body_force_dir[0]*node.X + body_force_dir[1]*node.Y + body_force_dir[2]*node.Z    # Iteration node body force projection
            rho = node.GetSolutionStepValue(KratosMultiphysics.DENSITY, 0)                                      # Nodal density value
            hyd_pres = rho*body_force_norm*(self.h_top + (body_force_proj-min_proj))                            # Iteration node hydrostatic pressure
            cur_pres = node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE, 0)                       # Iteration node imposed external pressure

            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, hyd_pres+cur_pres)                        # Add the hydrostatic component to the PRESSURE value
            node.SetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE, 0, hyd_pres+cur_pres)               # Add the hydrostatic component to the EXTERNAL_PRESSURE value


    def _ComputeOutletCharacteristicVelocity(self):
        # Compute the outlet average velocity
        outlet_avg_vel_norm = 0.0
        for node in self.outlet_model_part.GetCommunicator().LocalMesh().Nodes:
            vnode = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            outlet_avg_vel_norm += math.sqrt(vnode[0]*vnode[0] + vnode[1]*vnode[1] + vnode[2]*vnode[2])
        outlet_avg_vel_norm = self.outlet_model_part.GetCommunicator().SumAll(outlet_avg_vel_norm)

        tot_len = len(self.outlet_model_part.GetCommunicator().LocalMesh().Nodes)   # Partition outlet model part number of nodes
        tot_len = self.outlet_model_part.GetCommunicator().SumAll(tot_len)          # Get the total outlet model part nodes

        outlet_avg_vel_norm /= tot_len;

        # Store the average velocity in the ProcessInfo to be used in the outlet inflow prevention condition
        min_outlet_avg_vel_norm = 1.0
        if (outlet_avg_vel_norm >= min_outlet_avg_vel_norm):
            self.outlet_model_part.ProcessInfo[KratosFluid.CHARACTERISTIC_VELOCITY] = outlet_avg_vel_norm
        else:
            self.outlet_model_part.ProcessInfo[KratosFluid.CHARACTERISTIC_VELOCITY] = min_outlet_avg_vel_norm
