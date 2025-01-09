from math import sqrt
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsHydraulicsApplication as KratosFluidHydraulics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluidDynamics
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility
import KratosMultiphysics.kratos_utilities as KratosUtils

# def Factory(settings, Model):
#     if(type(settings) != KratosMultiphysics.Parameters):
#         raise Exception("expected input shall be a Parameters object, encapsulating a json string")
#     return ApplyHydraulicInletProcess(Model, settings["Parameters"])


class ApplyHydraulicInletProcess(KratosMultiphysics.Process):
    """
    Class to compute hydraulic inlet profile
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
        self.inlet_model_part = Model[settings["inlet_model_part_name"].GetString()]
        variable_name =settings["water_depth_variable"].GetString()
        self.water_depth_variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
        self.critical_depth_tolerance = settings["critical_depth_tolerance"].GetDouble()
        self.water_depth_tolerance = settings["water_depth_tolerance"].GetDouble()
        self.maximum_iterations = settings["maximum_iterations"].GetDouble()
        self.gravity = settings["gravity"].GetDouble()

        # Set the input vale type data
        self.discharge_value_is_constant= False
        self.discharge_value_is_function = False
        self.discharge_value_is_table = False

        if settings["value"].IsNumber():
            self.discharge_value_is_constant = True
            self.inlet_discharge = settings["value"].GetDouble()
        elif settings["value"].IsString():
            self.discharge_value_is_function = True
            self.function_string = settings["value"].GetString()
            self.aux_function = KratosMultiphysics.GenericFunctionUtility(self.function_string)
        else:
            self.discharge_value_is_table = True

            self.table = ReadCsvTableUtility(settings["value"]).Read(self.inlet_model_part)

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
            "value":{},
            "interval"        : [0.0,"End"],
            "water_depth_variable": "AUX_DISTANCE",
            "water_depth_tolerance":1e-8,
            "critical_depth_tolerance":1e-8,
            "maximum_iterations":100,
            "gravity"      : 9.81
                                    
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

        #Calculate the initial water depth guess.
        if not self.initial_water_depth:
            self.initial_water_depth = KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.InitialWaterDepth(
                self.inlet_model_part)

        # Determine the critical water depth based on the given inlet discharge value and the shape of the inlet model part. It is an iterative process based on the bisection method.

        # Define an interval within which the critical water depth belongs.
        # TODO: There is a more optimal way to define this interval in order to avoid iterations.
        water_depth_a = 1e-6
        water_depth_b = max([node.Z for node in self.inlet_model_part.Nodes])
        # The froude number for supercritical(small water depth) flows has a big value while for subcritical(High water depth) flows the froude number is small less than one.
        froude_number_a = 100000
        froude_number_c = 0.0

        #  Bisection method.
        iter_count = 0
        while (abs(froude_number_c - 1) > self.critical_depth_tolerance) and iter_count < self.maximum_iterations:
            water_depth_c = 0.5*(water_depth_a+water_depth_b)
            froude_number_c = self.CalculateFroudeNumber(water_depth_c)
            if (froude_number_c-1)*(froude_number_a-1)<0:
                water_depth_b=water_depth_c
            else:
                water_depth_a=water_depth_c
                froude_number_a = froude_number_c
            iter_count +=1
        if iter_count == self.maximum_iterations:
            KratosMultiphysics.Logger.PrintWarning("::[KratosFluidHydraulics]::", "Maximum iteractions for calculating critical water depth have been exceeded.")

        # Obtain the corresponding velocity according to the wetted area and the inlet discharge value.
        inlet_velocity = self.inlet_discharge/self.wetted_area

        # Assing the inlet velocity to inlet nodes and fix it.
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplexNonHistorical(self.inlet_model_part,self.domain_size,KratosFluidDynamics.INLET_NORMAL)
        KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.SetInletVelocity(
            self.inlet_model_part, inlet_velocity, self.water_depth_variable)

        #Assign the identical value of the inlet water depth to the DISTANCE variable (free surface) for all nodes associated with the inlet model part.
        KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.SetInletFreeSurface(
            self.inlet_model_part, KratosMultiphysics.INLET, self.water_depth_variable)


    def ExecuteFinalizeSolutionStep(self):
        # Here we free all of the nodes in the inlet
        if self.step_is_active:
            KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.FreeInlet(
                self.inlet_model_part)
        self.step_is_active = False

    def CalculateFroudeNumber(self, water_depth):
        #The Froude number is computed, and this value will be included into the objective function for the bisection method.
        # Assing the inlet free surface
        for node in self.inlet_model_part.Nodes:
            aux_distance = node.Z-water_depth
            if abs(aux_distance) < self.water_depth_tolerance:
                aux_distance = self.water_depth_tolerance if aux_distance > 0.0 else - self.water_depth_tolerance
            node.SetValue(self.water_depth_variable, aux_distance)

        # Calculate the wetted area and perimeter
        self.wetted_area = KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.CalculateWettedArea(
            self.inlet_model_part, KratosMultiphysics.INLET, self.water_depth_variable, False)
        self.wetted_perimeter = KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.CalculateWettedPetimeter(
            self.inlet_model_part, KratosMultiphysics.INLET, self.water_depth_variable, False)

        # Calculate froude number.
        froude_number = self.inlet_discharge/self.wetted_area /(sqrt(self.gravity * self.wetted_area / self.wetted_perimeter))
        return froude_number



