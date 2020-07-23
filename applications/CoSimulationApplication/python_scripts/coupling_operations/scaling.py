from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.function_callback_utility import GenericCallFunction

def Create(*args):
    return ScalingOperation(*args)

class ScalingOperation(CoSimulationCouplingOperation):
    """This operation performs scaling of values on an InterfaceData
    The value can be given directly as a value or as a string containing an evaluable function
    """
    def __init__(self, settings, solver_wrappers, process_info):
        if not settings.Has("scaling_factor"):
            raise Exception('Please provide a "scaling_factor"!')

        if settings["scaling_factor"].IsString():
            self.scaling_factor = settings["scaling_factor"].GetString()
        elif settings["scaling_factor"].IsNumber():
            self.scaling_factor = settings["scaling_factor"].GetDouble()
        else:
            raise Exception('The "scaling_factor" can only be provided as a number or a function-string')

        # removing since the type of "scaling_factor" can be double or string and hence would fail in the validation
        settings.RemoveValue("scaling_factor")

        super(ScalingOperation, self).__init__(settings, process_info)

        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)

    def Execute(self):
        process_info = self.interface_data.GetModelPart().ProcessInfo
        time = process_info[KM.TIME]
        step = process_info[KM.STEP]

        if not KM.IntervalUtility(self.settings).IsInInterval(time):
            if self.echo_level > 0:
                cs_tools.cs_print_info("ScalingOperation", "Skipped, not in interval")
            return

        if isinstance(self.scaling_factor, str):
            # TODO maybe make this to use COSIM_TIME and COSIM_STEP, such that it is generic for all solvers
            scope_vars = {'t' : time, 'step' : step} # make time and step useable in function
            current_scaling_factor = GenericCallFunction(self.scaling_factor, scope_vars, check=False) # evaluating function string
        else:
            current_scaling_factor = self.scaling_factor

        if self.echo_level > 0:
            cs_tools.cs_print_info("ScalingOperation", "Scaling-Factor", current_scaling_factor)
        self.interface_data.SetData(current_scaling_factor*self.interface_data.GetData()) # setting the scaled data

    def Check(self):
        if isinstance(self.scaling_factor, str):
            scope_vars = {'t' : 0.1, 'step' : 2} # dummy for testing
            GenericCallFunction(self.scaling_factor, scope_vars, check=True) # trying to evaluate the function string, such that the check can be disabled later

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED",
            "interval"  : [0.0, 1e30]
        }""")
        this_defaults.AddMissingParameters(super(ScalingOperation, cls)._GetDefaultSettings())
        return this_defaults
