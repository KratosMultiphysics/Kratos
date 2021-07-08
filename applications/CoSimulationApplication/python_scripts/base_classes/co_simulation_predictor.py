# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# other imports
from collections import deque

class CoSimulationPredictor:
    """Baseclass for the predictors used for CoSimulation
    It predicts the solution of the next step at the beginning of a step
    """
    def __init__(self, settings, solver_wrapper):
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(self._GetDefaultParameters())

        self.interface_data = solver_wrapper.GetInterfaceData(self.settings["data_name"].GetString())
        self.historical_data_accessor = HistoricalDataAccessor(self.interface_data, self._GetMinimumBufferSize())

        self.echo_level = self.settings["echo_level"].GetInt()

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        self.historical_data_accessor.CloneTimeStep()

    def Predict(self):
        raise Exception('"Predict" has to be implemented in the derived class!')

    def FinalizeSolutionStep(self):
        pass

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        cs_tools.cs_print_info("Predictor", colors.bold(self._ClassName()))

    def Check(self):
        print("The predictors do not yet implement Check!")

    def _UpdateData(self, updated_data):
        self.interface_data.SetData(updated_data)

        if self.echo_level > 3:
            cs_tools.cs_print_info(self._ClassName(), "Computed prediction")

    # returns the buffer size needed by the predictor. Can be overridden in derived classes
    def _GetMinimumBufferSize(self):
        return 2

    @classmethod
    def _ClassName(cls):
        return cls.__name__

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "type"       : "UNSPECIFIED",
            "solver"     : "UNSPECIFIED",
            "data_name"  : "UNSPECIFIED",
            "echo_level" : 0
        }""")


class HistoricalDataAccessor:
    def __init__(self, interface_data, buffer_size):
        self.interface_data = interface_data
        self.buffer_size = buffer_size
        self.step = 0

        self.additional_buffer_size = max(0, buffer_size - interface_data.GetBufferSize())
        self.aux_data = deque(maxlen=self.additional_buffer_size)

    def CloneTimeStep(self):
        self.step += 1
        if self.additional_buffer_size > 0:
            self.aux_data.appendleft(self.interface_data.GetData())

    def GetData(self, buffer_index):
        if buffer_index >= self.buffer_size:
            raise Exception("the buffer-size is not large enough current buffer size: {} | requested solution_step_index: {}!".format(self.buffer_size, buffer_index+1))

        if self.interface_data.GetBufferSize() > buffer_index:
            # the interface data buffer size is large enough
            return self.interface_data.GetData(buffer_index)
        else:
            # TODO check if time buffer is initialized
            return self.aux_data[buffer_index-self.additional_buffer_size]

    def TimeBufferIsInitialized(self):
        return self.step >= self.additional_buffer_size
