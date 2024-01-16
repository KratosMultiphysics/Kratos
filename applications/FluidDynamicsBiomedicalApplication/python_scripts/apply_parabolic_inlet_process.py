# Importing the Kratos Library
from typing import Type
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility

# Import applications
import KratosMultiphysics.FluidDynamicsBiomedicalApplication as KratosBio

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyParabolicInletProcess(model, settings["Parameters"])

class ApplyParabolicInletProcess(KratosMultiphysics.Process):
    """
    Class to compute parabolic inlet profile
    """

    def __init__(self, model, settings):
        # Base class constructor call
        KratosMultiphysics.Process.__init__(self)

        # Set default settings
        # Note the trick to allow scalar, function and table values of the parabola maximum value
        default_settings = self.GetDefaultParameters()
        if settings.Has("value"):
            if settings["value"].IsString():
                default_settings["value"].SetString("0.0")
            elif settings["value"].IsDouble():
                default_settings["value"].SetDouble(0.0)
        else:
            raise Exception("'value' not found. It needs to be user-provided.")

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # Check default settings
        settings.ValidateAndAssignDefaults(default_settings)

        # Check user-provided data
        if not settings["inlet_model_part_name"].GetString():
            raise ValueError("'inlet_model_part' not provided.")

        if not settings["wall_model_part_name"].GetString():
            raise ValueError("'wall_model_part' needs to be provided for 'parabolic' inlet distribution.")

        # Set the maximum value input data
        self.max_value_is_numeric = False
        self.max_value_is_function = False
        self.max_value_is_table = False
        if settings["value"].IsNumber():
            self.max_value_is_numeric = True
            self.max_value = settings["value"].GetDouble()
        elif settings["value"].IsString():
            self.max_value_is_function = True
            self.function_string = settings["value"].GetString()
            self.max_value_function = KratosMultiphysics.GenericFunctionUtility(self.function_string, settings["local_axes"])
        else:
            self.max_value_is_table = True
            inlet_model_part = model.GetModelPart(settings["inlet_model_part_name"].GetString())
            self.table = ReadCsvTableUtility(settings["value"]).Read(inlet_model_part)
        self.value_is_average = settings["value_is_average"].GetBool()
        self.value_is_flow_rate = settings["value_is_flow_rate"].GetBool()

        # Save model and settings containers
        self.model = model
        self.settings = settings

    @staticmethod
    def GetDefaultParameters():
        default_settings = KratosMultiphysics.Parameters("""{
            "wall_model_part_name": "",
            "inlet_model_part_name": "",
            "value" : {},
            "value_is_average" : true,
            "value_is_flow_rate" : true,
            "interval" : [0.0,"End"],
            "local_axes" : {},
            "parallel_distance_max_levels" : 25
        }""")
        return default_settings

    def ExecuteInitialize(self):
        # Get and check domain size
        inlet_model_part = self.model.GetModelPart(self.settings["inlet_model_part_name"].GetString())
        domain_size = inlet_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2,3]:
            raise ValueError(f"Wrong 'DOMAIN_SIZE' value {domain_size} in ProcessInfo container.")

        # Set the INLET flag in the inlet model part nodes and conditions
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INLET, True, inlet_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INLET, True, inlet_model_part.Conditions)

        # Compute the normal on the nodes of interest
        # Note that a custom normal variable is use to avoid interfering with the wall one
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplexNonHistorical(
            inlet_model_part,
            domain_size,
            KratosCFD.INLET_NORMAL)

        # Create an auxiliary volumetric model part with the elements attached to the inlet
        # On top of accelerating the wall distance calculation, this prevents missbehaviors in presence of complex geometries
        aux_inlet_model_part = KratosBio.ParabolicProfileUtilities.CreateAndFillInletAuxiliaryVolumeModelPart(inlet_model_part)

        # Prepare skin for wall distance calculation
        max_levels = self.settings["parallel_distance_max_levels"].GetInt()
        wall_model_part = self.model.GetModelPart(self.settings["wall_model_part_name"].GetString())
        KratosBio.ParabolicProfileUtilities.CalculateWallParallelDistance(wall_model_part, aux_inlet_model_part, max_levels)

        # Calculate the inlet area to do the flow rate to velocity conversion
        if self.value_is_flow_rate:
            inlet_model_part = self.model.GetModelPart(self.settings["inlet_model_part_name"].GetString())
            self.inlet_area = KratosBio.ParabolicProfileUtilities.CalculateInletArea(inlet_model_part)
            if self.inlet_area < 1.0e-12:
                KratosMultiphysics.Logger.PrintWarning(f"Inlet area {self.inlet_area} is close to zero in model part {inlet_model_part.FullName()}")

    def ExecuteBeforeSolutionLoop(self):
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        # Set the max value factor as the area quotient if the max value is a flow rate
        # Otherwise, if the max value is already a velocity set the max value factor to one
        # Also check if the provided value is the peak or the average value of the parabola
        value_is_average_factor = 2.0 if self.value_is_average else 1.0
        if self.value_is_flow_rate:
            max_value_factor = value_is_average_factor / self.inlet_area
        else:
            max_value_factor = value_is_average_factor

        # Set the parabolic inlet values
        inlet_model_part = self.model.GetModelPart(self.settings["inlet_model_part_name"].GetString())
        current_time = inlet_model_part.ProcessInfo[KratosMultiphysics.TIME]
        if self.interval.IsInInterval(current_time):
            self.step_is_active = True
            if self.max_value_is_numeric:
                KratosBio.ParabolicProfileUtilities.ImposeParabolicInlet(inlet_model_part, self.max_value, max_value_factor)
            elif self.max_value_is_function:
                KratosBio.ParabolicProfileUtilities.ImposeParabolicInlet(inlet_model_part, self.max_value_function, max_value_factor)
            elif self.max_value_is_table:
                current_max_value = self.table.GetValue(current_time)
                KratosBio.ParabolicProfileUtilities.ImposeParabolicInlet(inlet_model_part, current_max_value, max_value_factor)
            else:
                raise TypeError("Wrong maximum value data type.")

    def ExecuteFinalizeSolutionStep(self):
        # Here we free all of the nodes in the inlet
        if self.step_is_active:
            inlet_model_part = self.model.GetModelPart(self.settings["inlet_model_part_name"].GetString())
            KratosBio.ParabolicProfileUtilities.FreeParabolicInlet(inlet_model_part)
        self.step_is_active = False

    @classmethod
    def _ReturnDistanceProcessPrototype(cls, domain_size):
        if domain_size == 2:
            return KratosMultiphysics.ParallelDistanceCalculationProcess2D
        else:
            return KratosMultiphysics.ParallelDistanceCalculationProcess3D
