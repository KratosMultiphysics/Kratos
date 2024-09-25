import math
import KratosMultiphysics
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

# Import applications
import KratosMultiphysics.FluidDynamicsBiomedicalApplication as KratosBio


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyWindkesselOutletProcess(Model, settings["Parameters"])


class ApplyWindkesselOutletProcess(KratosMultiphysics.Process):
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
            "resistance1"        : 0.0,
            "resistance2"        : 0.0,
            "compliance"         : 0.0,
            "initial_pressure"   : 0.0
        }
        """)

        # Trick: allows "value" to be a double, a string or a table value (otherwise the ValidateAndAssignDefaults might fail)
        if(settings.Has("value")):
            if(settings["value"].IsString()):
                default_settings["value"].SetString("0.0")
            elif settings["value"].IsNumber():
                default_settings["value"].SetDouble(0.0)
        else:
            err_msg = "Provided settings have no 'value'. This needs to be provided."
            raise Exception(err_msg)

        settings.ValidateAndAssignDefaults(default_settings)

        # Set a Kratos parameters suitable for the core processes to set the PRESSURE
        pres_settings = settings.Clone()
        pres_settings.RemoveValue("resistance1")
        pres_settings.RemoveValue("resistance2")
        pres_settings.RemoveValue("compliance")
        pres_settings.RemoveValue("initial_pressure")

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

        self.R1 = settings["resistance1"].GetDouble()
        self.R2 = settings["resistance2"].GetDouble()
        self.C  = settings["compliance"].GetDouble()
        p0_mmHg = settings["initial_pressure"].GetDouble()

        self.conv = 13.545*9.81
        p0 = p0_mmHg*self.conv

        self.previous_q1 = 0.0
        self.current_p1 = p0

        # Set the OUTLET flag in the outlet model part nodes and conditions
        self.outlet_model_part = Model[pres_settings["model_part_name"].GetString()]
        for node in self.outlet_model_part.Nodes:
            node.Set(KratosMultiphysics.OUTLET, True)
        for condition in self.outlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.OUTLET, True)


        #self.pres_settings = pres_settings
        #self.ext_pres_settings = ext_pres_settings

        # Construct the base process AssignValueProcess
        #self.aux_pressure_process = AssignScalarVariableProcess(Model, self.pres_settings)
        #self.aux_external_pressure_process = AssignScalarVariableProcess(Model, self.ext_pres_settings)



    def ExecuteInitializeSolutionStep(self):

        # Here the value to be provided to the pressure is computed from an external calculation
        delta_t = self.outlet_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        t = self.outlet_model_part.ProcessInfo[KratosMultiphysics.TIME]

        current_q1 = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRate(self.outlet_model_part)
        print('Current flow rate', current_q1)

        modified_p1 = (1/self.C*(current_q1*( 1 + self.R1/self.R2) + self.R1*self.C*(current_q1 - self.previous_q1)/delta_t - self.current_p1/self.R2))*delta_t + self.current_p1
        print('Outlet new pressure:', modified_p1/self.conv)

        self.pres_settings["value"].SetDouble(modified_p1)
        self.ext_pres_settings["value"].SetDouble(modified_p1)

        # Call the base process ExecuteInitializeSolutionStep()
        self.aux_pressure_process.ExecuteInitializeSolutionStep()
        self.aux_external_pressure_process.ExecuteInitializeSolutionStep()

        self.previous_q1 = current_q1
        self.current_p1 = modified_p1


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
        min_proj = 1.0e15
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


    def __ComputeOutletAverageVelocity(self):
        # Compute the outlet average velocity
        outlet_avg_vel_norm = 0.0
        for node in self.outlet_model_part.GetCommunicator().LocalMesh().Nodes:
            vnode = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            outlet_avg_vel_norm += math.sqrt(vnode[0]*vnode[0] + vnode[1]*vnode[1] + vnode[2]*vnode[2])
        comm = self.outlet_model_part.GetCommunicator().GetDataCommunicator()
        outlet_avg_vel_norm = comm.SumAll(outlet_avg_vel_norm)

        tot_len = len(self.outlet_model_part.GetCommunicator().LocalMesh().Nodes)   # Partition outlet model part number of nodes
        tot_len = comm.SumAll(tot_len)                                              # Get the total outlet model part nodes

        outlet_avg_vel_norm /= tot_len

        # Check the outlet average velocity minimum value and return
        min_outlet_avg_vel_norm = 1.0e-12
        return outlet_avg_vel_norm if outlet_avg_vel_norm >= min_outlet_avg_vel_norm else min_outlet_avg_vel_norm

    def __ComputeOutletMaxVelocity(self):
        # Compute the outlet max velocity
        outlet_max_vel_norm = 0.0
        for node in self.outlet_model_part.GetCommunicator().LocalMesh().Nodes:
            v_node = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            aux_val = v_node[0]*v_node[0] + v_node[1]*v_node[1] + v_node[2]*v_node[2]
            if aux_val > outlet_max_vel_norm:
                outlet_max_vel_norm = aux_val
        outlet_max_vel_norm = self.outlet_model_part.GetCommunicator().GetDataCommunicator().MaxAll(outlet_max_vel_norm)
        outlet_max_vel_norm = math.sqrt(outlet_max_vel_norm)

        # Check the outlet max velocity minimum value and return
        min_outlet_vel_norm = 1.0e-12
        return outlet_max_vel_norm if outlet_max_vel_norm >= min_outlet_vel_norm else min_outlet_vel_norm
