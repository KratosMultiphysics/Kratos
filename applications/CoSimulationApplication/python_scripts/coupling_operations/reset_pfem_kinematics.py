# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# Additional imports

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(*args):
    return ResetPfemKinematics(*args)

class ResetPfemKinematics(CoSimulationCouplingOperation):
    """This operation is used to reset the PFEM kinematic values on a given model part.
    It does not touch the SOLID/RIGID nodes
    (See the PfemFluidDynamicsApp - move_mesh_utility.cpp implementation)
    TODO:
    - add tests
    - more cleanup
    """
    def __init__(self, settings, solver_wrappers, process_info):
        super().__init__(settings, process_info)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.model_part = self.model[self.model_part_name]

        # Detect 'End' as a tag and replace it by a large number
        if(self.settings.Has('interval')):
            if(self.settings['interval'][1].IsString()):
                if(self.settings['interval'][1].GetString() == 'End' or self.settings['interval'][1].GetString() == 'end'):
                    self.settings['interval'][1].SetDouble(1e30)
                else:
                    raise Exception('The second value of interval can be \'End\' or a number, interval currently:' + self.settings['interval'].PrettyPrintJsonString())
        self.interval = self.settings["interval"].GetVector()

    def InitializeCouplingIteration(self):
        current_time = self.model_part.ProcessInfo[KM.TIME]
        if((current_time >= self.interval[0]) and (current_time < self.interval[1])):
            self._ResetPfemKinematicValues()

    def _ResetPfemKinematicValues(self):
        KM.PfemFluidDynamicsApplication.MoveMeshUtility().ResetPfemKinematicValues(self.model_part)

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                : "UNSPECIFIED",
            "model_part_name"       : "",
            "interval"              : [0.0, 1e30]
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
