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
            "initial_pressure"   : 0.0,
            "pressure_in_mmHg"   : true,
            "echo_level"         : 0
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
        pres_settings.RemoveValue("echo_level")

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

        self.R1      = settings["resistance1"].GetDouble()
        self.R2      = settings["resistance2"].GetDouble()
        self.C       = settings["compliance"].GetDouble()
        p0_mmHg      = settings["initial_pressure"].GetDouble()
        self.echo    = settings["echo_level"].GetInt()
        self.pressure_in_mmHg = settings["pressure_in_mmHg"].GetBool()

        self.conv = 13.545*9.81 # Pressure conversion factor. It is used if pressure is provided in mmHg
        if self.pressure_in_mmHg is True:
            p0 = p0_mmHg*self.conv
        else:
            p0 = p0_mmHg   # The pressure variable passed from the json file is given in Pa


        self.previous_q1 = 0.0
        self.current_p1 = p0 # in Pa

        # Set the OUTLET flag in the outlet model part nodes and conditions
        self.outlet_model_part = Model[pres_settings["model_part_name"].GetString()]
        for node in self.outlet_model_part.Nodes:
            node.Set(KratosMultiphysics.OUTLET, True)
        for condition in self.outlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.OUTLET, True)

    def ExecuteInitializeSolutionStep(self):

        # Here the value to be provided to the outlet pressure is computed as the result of an ODE:
        delta_t = self.outlet_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        t = self.outlet_model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.current_q1 = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRate(self.outlet_model_part)

        self.modified_p1 = (1/self.C*(self.current_q1*( 1 + self.R1/self.R2) + self.R1*self.C*(self.current_q1 - self.previous_q1)/delta_t - self.current_p1/self.R2))*delta_t + self.current_p1

        if self.echo > 0:
            print('Current flow rate', self.current_q1)
            if self.pressure_in_mmHg:
                print('Outlet new pressure:', self.modified_p1/self.conv, " mmHg")
            else:
                print('Outlet new pressure:', self.modified_p1, " Pa")

        for node in self.outlet_model_part.Nodes:
            # Setting new solution on the nodes
            node.Fix(KratosMultiphysics.PRESSURE)
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,self.modified_p1)
            node.SetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,self.modified_p1)

    def ExecuteFinalizeSolutionStep(self):

        self.previous_q1 = self.current_q1
        self.current_p1 = self.modified_p1


    # Private methods section
