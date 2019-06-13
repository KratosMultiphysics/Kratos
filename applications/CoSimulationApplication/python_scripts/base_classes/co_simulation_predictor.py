from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class CosimulationPredictor(object):
    def __init__(self, settings, solver):
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(self._GetDefaultSettings())

        self.solver = solver
        self.interface_data = self.solver.GetInterfaceData(self.settings["data_name"].GetString())

        self.echo_level = self.settings["echo_level"].GetInt()

        # TODO check buffer size
        self._GetMinimumBufferSize()

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def Predict(self):
        raise Exception('"Predict" has to be implemented in the derived class!')

    def FinalizeSolutionStep(self):
        pass

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        cs_tools.classprint("Predictor", cs_tools.bold(self._Name()))

    def Check(self):
        print("The predictors do not yet implement Check!")

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        return self.__class__.__name__

    def _UpdateData(self, updated_data):
        self.interface_data.ApplyUpdateToData(updated_data)

        if self.echo_level > 3:
            cs_tools.classprint(self._Name(), "Computed prediction")


    # returns the buffer size needed by the predictor. Can be overridden in derived classes
    def _GetMinimumBufferSize(self):
        return 2

    @classmethod
    def _GetDefaultSettings(cls):
        return cs_tools.cs_data_structure.Parameters("""{
            "type"       : "UNSPECIFIED",
            "solver"     : "UNSPECIFIED",
            "data_name"  : "UNSPECIFIED",
            "echo_level" : 0
        }""")
