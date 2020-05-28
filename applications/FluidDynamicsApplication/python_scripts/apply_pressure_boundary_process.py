import math
import KratosMultiphysics
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyPressureBoundaryProcess(Model, settings["Parameters"])


class ApplyPressureBoundaryProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "mesh_id"            : 0,
            "model_part_name"    : "",
            "variable_name"      : "PRESSURE",
            "value"              : 0.0,
            "constrained"        : true,
            "interval"           : [0.0,"End"],
            "hydrostatic"        : false,
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
        pres_settings.RemoveValue("hydrostatic")
        pres_settings.RemoveValue("h_top")

        # Create a copy of the PRESSURE settings to set the EXTERNAL_PRESSURE
        ext_pres_settings = pres_settings.Clone()
        ext_pres_settings["constrained"].SetBool(False)
        ext_pres_settings["variable_name"].SetString("EXTERNAL_PRESSURE")

        # Check the core processes input data
        if (pres_settings["model_part_name"].GetString() == ""):
            raise Exception("Empty pressure boundary model part name. Set a valid model part name.")
        elif (pres_settings["value"].IsString()):
            if (pres_settings["value"].GetString == ""):
                raise Exception("Boundary pressure value is empty.")

        self.hydrostatic = settings["hydrostatic"].GetBool()
        self.h_top = settings["h_top"].GetDouble()

        # Set the PRESSURE_BC flag in the outlet model part nodes and conditions
        self.outlet_model_part = Model[pres_settings["model_part_name"].GetString()]
        for node in self.outlet_model_part.Nodes:
            node.Set(KratosMultiphysics.PRESSURE_BC, True)
        for condition in self.outlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.PRESSURE_BC, True)

        # Construct the base process AssignValueProcess
        self.aux_pressure_process = AssignScalarVariableProcess(Model, pres_settings)
        self.aux_external_pressure_process = AssignScalarVariableProcess(Model, ext_pres_settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        self.aux_pressure_process.ExecuteInitializeSolutionStep()
        self.aux_external_pressure_process.ExecuteInitializeSolutionStep()

        # If considered, add the hydrostatic component to the boundary pressure
        if (self.hydrostatic):
            self._AddHydrostaticComponent()


    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        self.aux_pressure_process.ExecuteFinalizeSolutionStep()
        self.aux_external_pressure_process.ExecuteFinalizeSolutionStep()


    # Private methods section
    def _AddHydrostaticComponent(self):
        # Initialize body force value (avoid segfault in MPI if the local mesh has no pressure boundary nodes)
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
