import KratosMultiphysics

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
            "model_part_name"    : "",
            "variable_name"      : "PRESSURE",
            "constrained"        : true,
            "value"              : 0.0,
            "interval"           : [0.0,"End"],
            "characteristic_resistance"        : 0.0,
            "peripheral_resistance"        : 0.0,
            "arterial_compliance"         : 0.0,
            "venous_pressure"    : 0.0,
            "initial_pressure"   : 0.0,
            "pressure_unit"      : "mmHg",
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

        settings.ValidateAndAssignDefaults(settings)

        # Check the core processes input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty outlet pressure model part name. Set a valid model part name.")
        elif (settings["variable_name"].GetString() != "PRESSURE"):
            raise Exception("Outlet pressure settings variable_name is not PRESSURE.")
        elif (settings["value"].IsString()):
            if (settings["value"].GetString == ""):
                raise Exception("Outlet pressure function sting is empty.")
                raise Exception("Outlet external pressure function sting is empty.")

        self.R1      = settings["characteristic_resistance"].GetDouble()
        self.R2      = settings["peripheral_resistance"].GetDouble()
        self.C       = settings["arterial_compliance"].GetDouble()
        p0_mmHg      = settings["initial_pressure"].GetDouble()
        pv_mmHg      = settings["venous_pressure"].GetDouble()
        self.echo    = settings["echo_level"].GetInt()
        self.pressure_unit = settings["pressure_unit"].GetString()

        self.conv = 13.545*9.81 # Pressure conversion factor. It is used if pressure is provided in mmHg
        if self.pressure_unit == "mmHg" :
            self.pv = pv_mmHg*self.conv
            p0 = p0_mmHg*self.conv
        elif self.pressure_unit == "Pa" :
            self.pv = pv_mmHg
            p0      = p0_mmHg   # The pressure variable passed from the json file is given in Pa
        else :
            raise Exception("Pressure unit measyre can be given in mmHg or in Pa")

        self.previous_q1 = 0.0
        self.current_p1 = p0 # in Pa

        # Set the OUTLET flag in the outlet model part nodes and conditions
        self.outlet_model_part = Model[settings["model_part_name"].GetString()]
        for node in self.outlet_model_part.Nodes:
            node.Set(KratosMultiphysics.OUTLET, True)
        for condition in self.outlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.OUTLET, True)

    def ExecuteInitializeSolutionStep(self):

        # Here the value to be provided to the outlet pressure is computed as the result of an ODE:
        delta_t = self.outlet_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        t = self.outlet_model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.current_q1 = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRate(self.outlet_model_part)

        self.modified_p1 = (1/self.C*(self.current_q1*( 1 + self.R1/self.R2) + self.R1*self.C*(self.current_q1 - self.previous_q1)/delta_t - (self.current_p1 - self.pv)/self.R2))*delta_t + self.current_p1

        if self.echo > 0:
            KratosMultiphysics.Logger.PrintInfo("Windkessel", f"Current flow rate: {self.current_q1}")
            if self.pressure_unit == "mmHg" :
                KratosMultiphysics.Logger.PrintInfo("Windkessel", f"Outlet new pressure: {self.modified_p1/self.conv} mmHg")
            elif self.pressure_unit == "Pa" :
                KratosMultiphysics.Logger.PrintInfo("Windkessel", f"Outlet new pressure: {self.modified_p1} Pa")
            else :
                raise Exception("Pressure unit measure can be given in mmHg or in Pa")

        for node in self.outlet_model_part.Nodes:
            # Setting new solution on the nodes
            node.Fix(KratosMultiphysics.PRESSURE)
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,self.modified_p1)
            node.SetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,self.modified_p1)

    def ExecuteFinalizeSolutionStep(self):

        self.previous_q1 = self.current_q1
        self.current_p1 = self.modified_p1


    # Private methods section
