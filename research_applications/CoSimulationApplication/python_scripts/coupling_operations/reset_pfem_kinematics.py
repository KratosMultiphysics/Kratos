# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

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
    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        model = solver_wrappers[self.settings["solver"].GetString()].model
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.model_part = model[self.model_part_name]

        self.interval = KM.IntervalUtility(settings)

    def InitializeCouplingIteration(self):
        if self.interval.IsInInterval(self.model_part.ProcessInfo[KM.TIME]):
            self._ResetPfemKinematicValues()

            if self.echo_level > 0:
                cs_tools.cs_print_info(self._ClassName(), "PFEM KINEMATICS RESET IN MODEL PART: " + self.model_part_name)

    def _ResetPfemKinematicValues(self):
        KM.PfemFluidDynamicsApplication.PFEMMoveMeshUtility.ResetPfemKinematicValues(self.model_part)

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                : "UNSPECIFIED",
            "model_part_name"       : "",
            "interval"              : [0.0, 1e30]
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
