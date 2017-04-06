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
            "mesh_id"         : 0,
            "model_part_name" : "",
            "variable_name"   : "PRESSURE",
            "constrained"     : true,
            "value"           : 0.0
        }
        """)

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if(settings.Has("value")):
            if(settings["value"].IsString()):
                default_settings["value"].SetString("0.0")

        settings.ValidateAndAssignDefaults(default_settings)

        # Create a copy of the PRESSURE settings to set the EXTERNAL_PRESSURE
        ext_pres_settings = settings.Clone()
        ext_pres_settings["constrained"].SetBool(False)
        ext_pres_settings["variable_name"].SetString("EXTERNAL_PRESSURE")

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty outlet pressure model part name. Set a valid model part name.")
        elif (ext_pres_settings["model_part_name"].GetString() == ""):
            raise Exception("Empty outlet external pressure model part name. Set a valid model part name.")
        elif (settings["variable_name"].GetString() != "PRESSURE"):
            raise Exception("Outlet pressure settings variable_name is not PRESSURE.")
        elif (ext_pres_settings["variable_name"].GetString() != "EXTERNAL_PRESSURE"):
            raise Exception("Outlet external pressure settings variable_name is not EXTERNAL_PRESSURE.")
        elif (settings["value"].IsString()):
            if (settings["value"].GetString == ""):
                raise Exception("Outlet pressure function sting is empty.")
        elif (ext_pres_settings["value"].IsString()):
            if (ext_pres_settings["value"].GetString == ""):
                raise Exception("Outlet external pressure function sting is empty.")

        # Set the OUTLET flag in the outlet model part nodes and conditions
        self.outlet_model_part = Model[settings["model_part_name"].GetString()]
        for node in self.outlet_model_part.Nodes:
            node.Set(KratosMultiphysics.OUTLET, True)
        for condition in self.outlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.OUTLET, True)

        # Construct the base process AssignValueProcess
        import experimental_assign_value_process
        self.aux_pressure_process = experimental_assign_value_process.AssignValueProcess(Model, settings)
        self.aux_external_pressure_process = experimental_assign_value_process.AssignValueProcess(Model, ext_pres_settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        self.aux_pressure_process.ExecuteInitializeSolutionStep()
        self.aux_external_pressure_process.ExecuteInitializeSolutionStep()

        # Compute the current process average velocity
        outlet_avg_vel_norm = 0.0
        for node in self.outlet_model_part.GetCommunicator().LocalMesh().Nodes:
            vnode = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            outlet_avg_vel_norm += math.sqrt(vnode[0]**2+vnode[1]**2+vnode[2]**2)
        outlet_avg_vel_norm /= len(self.outlet_model_part.GetCommunicator().LocalMesh().Nodes)

        # Perform MPI communication
        self.outlet_model_part.GetCommunicator().SumAll(outlet_avg_vel_norm)                 # Accumulate the average velocity in each processor
        outlet_avg_vel_norm /= self.outlet_model_part.GetCommunicator().TotalProcesses()     # Compute the average between processors

        # Store the average velocity in the ProcessInfo to be used in the outlet inflow prevention condition
        min_outlet_avg_vel_norm = 1.0
        if (outlet_avg_vel_norm >= min_outlet_avg_vel_norm):
            self.outlet_model_part.GetRootModelPart().ProcessInfo[KratosFluid.CHARACTERISTIC_VELOCITY] = outlet_avg_vel_norm
        else:
            self.outlet_model_part.GetRootModelPart().ProcessInfo[KratosFluid.CHARACTERISTIC_VELOCITY] = min_outlet_avg_vel_norm


    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        self.aux_pressure_process.ExecuteFinalizeSolutionStep()
        self.aux_external_pressure_process.ExecuteFinalizeSolutionStep()
