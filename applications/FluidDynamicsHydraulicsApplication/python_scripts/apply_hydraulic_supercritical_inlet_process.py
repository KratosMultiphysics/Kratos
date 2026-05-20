from math import sqrt
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsHydraulicsApplication as KratosFluidHydraulics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluidDynamics
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility
import KratosMultiphysics.kratos_utilities as KratosUtils

# def Factory(settings, Model):
#     if(type(settings) != KratosMultiphysics.Parameters):
#         raise Exception("expected input shall be a Parameters object, encapsulating a json string")
#     return ApplyHydraulicSuperCriticalInletProcess(Model, settings["Parameters"])


class ApplyHydraulicSuperCriticalInletProcess(KratosMultiphysics.Process):
    """
    Class to compute hydraulic inlet profile
    """
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Set default settings
        # Note the trick to allow scalar, function and table values of the inlet discharge
        default_settings = self.GetDefaultParameters()
        if settings.Has("value_discharge"):
            if settings["value_discharge"].IsString():
                default_settings["value_discharge"].SetString("0.0")
            elif settings["value_discharge"].IsDouble():
                default_settings["value_discharge"].SetDouble(0.0)
        else:
            raise Exception("'value_discharge' not found. It needs to be user-provided.")

        if settings.Has("value_velocity"):
            if settings["value_velocity"].IsString():
                default_settings["value_velocity"].SetString("0.0")
            elif settings["value_velocity"].IsDouble():
                default_settings["value_velocity"].SetDouble(0.0)
        else:
            raise Exception(
                "'value_velocity' not found. It needs to be user-provided.")

        self.interval = KratosMultiphysics.IntervalUtility(settings)
        settings.ValidateAndAssignDefaults(default_settings)

        # Get the inlet model part , the inlet free surface (water depth) variable name and user tolerances.
        self.inlet_model_part = Model[settings["inlet_model_part_name"].GetString()]
        variable_name =settings["water_depth_variable"].GetString()
        self.tolerance = settings["tolerance"].GetString()
        self.maximum_iterations = settings["maximum_iterations"].GetString()
        self.water_depth_variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)

        # Set the input vale type data
        self.discharge_value_is_constant= False
        self.discharge_value_is_function = False
        self.discharge_value_is_table = False

        self.velocity_value_is_constant= False
        self.velocity_value_is_function = False
        self.velocity_value_is_table = False

        if settings["value_discharge"].IsNumber():
            self.discharge_value_is_constant = True
            self.inlet_discharge = settings["value_discharge"].GetDouble()
        elif settings["value_discharge"].IsString():
            self.discharge_value_is_function = True
            self.function_string = settings["value_discharge"].GetString()
            self.aux_function = KratosMultiphysics.GenericFunctionUtility(self.function_string)
        else:
            self.discharge_value_is_table = True

            self.table = ReadCsvTableUtility(
                settings["value_discharge"]).Read(self.inlet_model_part)

        if settings["value_velocity"].IsNumber():
            self.velocity_value_is_constant = True
            self.inlet_velocity = settings["value_velocity"].GetDouble()
        elif settings["value_discharge"].IsString():
            self.velocity_value_is_function = True
            self.function_string = settings["value_velocity"].GetString()
            self.aux_function = KratosMultiphysics.GenericFunctionUtility(
                self.function_string)
        else:
            self.velocity_value_is_table = True

            self.table = ReadCsvTableUtility(
                settings["value_velocity"]).Read(self.inlet_model_part)

        # Set INLET flag to all nodes and conditions belonging to inlet model part.
        for node in self.inlet_model_part.Nodes:
            node.Set(KratosMultiphysics.INLET, True)
        for condition in self.inlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.INLET, True)

        # Get and check domain size
        self.domain_size = self.inlet_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if self.domain_size not in [2,3]:
            raise ValueError(f"Wrong 'DOMAIN_SIZE' value {self.domain_size} in ProcessInfo container.")

        self.initial_water_depth = None

    @staticmethod
    def GetDefaultParameters():
        default_settings = KratosMultiphysics.Parameters("""
        {
            "inlet_model_part_name" : "",
            "value_discharge":{},
            "value_velocity":{},
            "interval"        : [0.0,"End"],
            "tolerance" :1e-3,
            "maximum_iterations" :10,
            "water_depth_variable": "AUX_DISTANCE"
        }
        """)
        return default_settings

    def ExecuteBeforeSolutionLoop(self):
        # Note that this is not required at all but we do it just in case any other class needs the inlet distance and velocity
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):

        # For each time step obtain the corresponding inlet discharge value according to the external data type; table, function or value.
        current_time = self.inlet_model_part.ProcessInfo[KratosMultiphysics.TIME]
        if self.interval.IsInInterval(current_time):
            self.step_is_active = True
            if self.discharge_value_is_function:

                if self.aux_function.DependsOnSpace() == False:  # depends on time only
                    self.inlet_discharge = self.aux_function.CallFunction(0.0, 0.0, 0.0, current_time, 0.0, 0.0, 0.0)
                else:  # Only time-dependent hydrograms are implemented.
                    err_msg = "There is no space-dependent function implemented "
                    raise ValueError(err_msg)

            elif self.discharge_value_is_table:
                self.inlet_discharge = self.table.GetValue(current_time)

        if self.interval.IsInInterval(current_time):
            self.step_is_active = True
            if self.velocity_value_is_function:

                if self.aux_function.DependsOnSpace() == False:  # depends on time only
                    self.inlet_velocity = self.aux_function.CallFunction(
                        0.0, 0.0, 0.0, current_time, 0.0, 0.0, 0.0)
                else:  # Only time-dependent hydrograms are implemented.
                    err_msg = "There is no space-dependent function implemented "
                    raise ValueError(err_msg)
            elif self.velocity_value_is_table:
                self.inlet_velocity = self.table.GetValue(current_time)


        #Calculate the initial water depth guess.
        if not self.initial_water_depth:
            self.initial_water_depth = KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.InitialWaterDepth(
                self.inlet_model_part)
        KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.SetBoundaryWaterDepth(self.inlet_model_part, self.initial_water_depth,  self.water_depth_variable)
        wetted_aux_area = KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.CalculateWettedArea(
            self.inlet_model_part, KratosMultiphysics.INLET, self.water_depth_variable, False)

        while abs(wetted_aux_area-self.real_wettted_area) > self.tolerance and iter_count < self.maximum_iterations:
            if self.wetted_area>self.real_wetted_area:
                self.initial_water_depth/=2
            else:
                2*self.real_wetted_area
            KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.SetBoundaryWaterDepth(
                    self.inlet_model_part, self.initial_water_depth,  self.water_depth_variable)
            wetted_aux_area = KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.CalculateWettedArea(
                self.inlet_model_part, KratosMultiphysics.INLET, self.water_depth_variable, False)
            iter_count += 1
            if iter_count == self.maximum_iterations:
                KratosMultiphysics.Logger.PrintWarning(
                "::[KratosFluidHydraulics]::", "Maximum iteractions for calculating wetted area have been exceeded.")



        # Assing the inlet velocity to inlet nodes and fix it.
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplexNonHistorical(self.inlet_model_part,self.domain_size,KratosFluidDynamics.INLET_NORMAL)
        KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.SetInletVelocity(
            self.inlet_model_part, self.inlet_velocity, self.water_depth_variable)

        #Assign the identical value of the inlet water depth to the DISTANCE variable (free surface) for all nodes associated with the inlet model part.
        KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.SetInletFreeSurface(
            self.inlet_model_part, KratosMultiphysics.INLET, self.water_depth_variable)


    def ExecuteFinalizeSolutionStep(self):
        # Here we free all of the nodes in the inlet
        if self.step_is_active:
            KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.FreeInlet(
                self.inlet_model_part)
        self.step_is_active = False



