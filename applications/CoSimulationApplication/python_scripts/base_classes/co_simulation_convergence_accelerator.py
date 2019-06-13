from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from . import co_simulation_coupling_operation
# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class CoSimulationConvergenceAccelerator(co_simulation_coupling_operation.CoSimulationCouplingOperation):
    def __init__(self, settings, solver):
        super(CoSimulationConvergenceAccelerator, self).__init__(settings)

        self.solver = solver
        self.interface_data = self.solver.GetInterfaceData(self.settings["data_name"].GetString())

    def InitializeCouplingIteration(self):
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update
        self.input_data = self.interface_data.GetNumpyArray()

    def ComputeAndApplyUpdate(self):
        current_data = self.interface_data.GetNumpyArray()
        residual = current_data - self.input_data
        updated_data = self.input_data + self.ComputeUpdate(residual, self.input_data)
        self.interface_data.ApplyUpdateToData(updated_data)

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        cs_tools.classprint("Convergence Accelerator", cs_tools.bold(self._Name()))

    def Check(self):
        print("ConvAcc does not yet implement Check")

    def Execute(self):
        raise Exception("This is not supposed to be used for ConvergenceAccelerators!")

    def ComputeUpdate( self, residual, previous_data ):
        raise Exception('"ComputeUpdate" has to be implemented in the derived class!')

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = cs_tools.cs_data_structure.Parameters("""{
            "type"       : "UNSPECIFIED",
            "solver"     : "UNSPECIFIED",
            "data_name"  : "UNSPECIFIED"
        }""")

        this_defaults.AddMissingParameters(super(CoSimulationConvergenceAccelerator, cls)._GetDefaultSettings())

        return this_defaults
