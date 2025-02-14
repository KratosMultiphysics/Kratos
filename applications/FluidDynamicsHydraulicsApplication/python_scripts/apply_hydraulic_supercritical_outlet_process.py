import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsHydraulicsApplication as KratosFluidHydraulics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluidDynamics
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility
import KratosMultiphysics.kratos_utilities as KratosUtils


class ApplyHydraulicSuperCriticalOutletProcess(KratosMultiphysics.Process):
    """
    Class to compute Supercritical Outlet Process
    """

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Set default settings
        # Note the trick to allow scalar, function and table values of the inlet discharge
        default_settings = self.GetDefaultParameters()
        if settings.Has("value"):
            if settings["value"].IsString():
                default_settings["value"].SetString("0.0")
            elif settings["value"].IsDouble():
                default_settings["value"].SetDouble(0.0)
        else:
            raise Exception("'value' not found. It needs to be user-provided.")

        self.interval = KratosMultiphysics.IntervalUtility(settings)
        settings.ValidateAndAssignDefaults(default_settings)

        # Get the inlet model part , the inlet free surface (water depth) variable name and user tolerances.
        self.outlet_model_part = Model[settings["outlet_model_part_name"].GetString(
        )]
        self.continuity_tolerance = settings["continuity_tolerance"].GetDouble(
        )
        self.inflow_detection = settings["inflow_detection"].GetBool()

        # Set the input vale type data
        self.discharge_value_is_constant = False
        self.discharge_value_is_function = False
        self.discharge_value_is_table = False

        if settings["value"].IsNumber():
            self.discharge_value_is_constant = True
            self.inlet_discharge = settings["value"].GetDouble()
        elif settings["value"].IsString():
            self.discharge_value_is_function = True
            self.function_string = settings["value"].GetString()
            self.aux_function = KratosMultiphysics.GenericFunctionUtility(
                self.function_string)
        else:
            self.discharge_value_is_table = True

            self.table = ReadCsvTableUtility(
                settings["value"]).Read(self.outlet_model_part)

        # Set OUTLET flag to all nodes and conditions belonging to inlet model part.
        for node in self.outlet_model_part.Nodes:
            node.Set(KratosMultiphysics.OUTLET, True)
        for condition in self.outlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.OUTLET, True)

        # Get and check domain size
        self.domain_size = self.outlet_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if self.domain_size not in [2, 3]:
            raise ValueError(
                f"Wrong 'DOMAIN_SIZE' value {self.domain_size} in ProcessInfo container.")

        self.initial_water_depth = None

    @staticmethod
    def GetDefaultParameters():
        default_settings = KratosMultiphysics.Parameters("""
        {
            "outlet_model_part_name" : "",
            "value":{},
            "interval"        : [0.0,"End"],
            "maximum_iterations":100,
            "continuity_tolerance":1e-1,
            "inflow_detection": true
        }
        """)
        return default_settings

    def ExecuteBeforeSolutionLoop(self):
        # Note that this is not required at all but we do it just in case any other class needs the inlet distance and velocity
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):

        # Obtain inlet dicharge for each time step
        current_time = self.outlet_model_part.ProcessInfo[KratosMultiphysics.TIME]
        if current_time < 0.0:
            self.outlet_model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.0
        if self.interval.IsInInterval(current_time):
            self.step_is_active = True
            if self.discharge_value_is_function:
                if self.aux_function.DependsOnSpace() == False:  # depends on time only
                    self.inlet_discharge = self.aux_function.CallFunction(
                        0.0, 0.0, 0.0, current_time, 0.0, 0.0, 0.0)
                else:  # Only time-dependent hydrograms are implemented.
                    err_msg = "There is no space-dependent function implemented "
                    raise ValueError(err_msg)

            elif self.discharge_value_is_table:
                self.inlet_discharge = self.table.GetValue(current_time)
        # Calculate Outlet water discharge value of the previous time step.
        outlet_water_discharge_n = KratosFluidDynamics.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(self.outlet_model_part)


        if abs(outlet_water_discharge_n - self.inlet_discharge) > self.continuity_tolerance:
            for node in self.outlet_model_part.Nodes:
                node.SetSolutionStepValue(
                    KratosMultiphysics.EXTERNAL_PRESSURE, 0.0)

        else:
            for node in self.outlet_model_part.Nodes:
                p_n = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE, 1)
                node.SetSolutionStepValue(
                    KratosMultiphysics.EXTERNAL_PRESSURE, p_n)

    def ExecuteFinalizeSolutionStep(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplexNonHistorical(
            self.outlet_model_part, self.domain_size, KratosFluidDynamics.OUTLET_NORMAL)
        if self.inflow_detection:
            KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.ApplyOutletInflowLimiter(
                self.outlet_model_part, KratosMultiphysics.VELOCITY, KratosFluidDynamics.OUTLET_NORMAL)
            # Is it necessary? or is it enough after solving vectorial convection problem?
            KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.ApplyOutletInflowLimiter(
                self.outlet_model_part, KratosFluidDynamics.FRACTIONAL_VELOCITY, KratosFluidDynamics.OUTLET_NORMAL)

