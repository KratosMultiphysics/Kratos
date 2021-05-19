# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

class CoSimulationPredictor:
    """Baseclass for the predictors used for CoSimulation
    It predicts the solution of the next step at the beginning of a step
    """
    def __init__(self, settings, solver_wrapper):
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(self._GetDefaultParameters())

        self.interface_data = solver_wrapper.GetInterfaceData(self.settings["data_name"].GetString())

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
