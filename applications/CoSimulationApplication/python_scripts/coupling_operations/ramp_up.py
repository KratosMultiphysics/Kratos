from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_wrappers):
    cs_tools.SettingsTypeCheck(settings)
    return RampUpOperation(settings, solver_wrappers)

class RampUpOperation(CoSimulationCouplingOperation):
    """This operation performs scaling of values on an InterfaceData
    The value can be given directly as a value or as a string containing an evaluable function
    """
    def __init__(self, settings, solver_wrappers):
        super(RampUpOperation, self).__init__(settings)

        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)

        self.num_steps = self.settings["num_steps"].GetInt()
        self.num_steps_executed = 0 # TODO this needs to be serialized otherwise will not work in restart

        ramp_type =  = self.settings["solver"].GetString()

        if ramp_type == "linear":
            self.ramp_function = lambda ....
        elif ramp_type == "sin^2":
            self.ramp_function = lambda ....
        else:
            raise Exception('The requested "ramp_type" "{}" is not available! Only the following options are available:\n\t"linear", "sin^2"'.format(ramp_type))

        #     return lambda conn: [conn[i] for i in [0, 2, 1, 3]]
        # elif salome_entity_type == "Hexa":
        #     return lambda conn: [conn[i] for i in [0, 3, 2, 1, 4, 7, 6, 5]]
        # elif salome_entity_type == "Penta":
        #     return lambda conn: [conn[i] for i in [0, 2, 1, 3, 5, 4]]
        # else:
        #     return lambda conn: conn

    def Execute(self):
        self.num_steps_executed += 1

        current_scaling_factor = self.ramp_function(self.num_steps_executed, self.num_steps)

        if self.echo_level > 0:
            cs_tools.cs_print_info("RampUpOperation", "Scaling-Factor", current_scaling_factor)
        # TODO skip if the data is restarted. This is an ugly fix and might fail if restart happesn earlier, but this is necessary until restart is implemented in CoSim
        self.interface_data.SetData(current_scaling_factor*self.interface_data.GetData()) # setting the scaled data

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED",
            "num_steps" : 10
            "ramp_type" : "sin^2"
        }""")
        this_defaults.AddMissingParameters(super(RampUpOperation, cls)._GetDefaultSettings())
        return this_defaults
