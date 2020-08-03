# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

def Create(*args):
    return ComputeNormalsOperation(*args)

class ComputeNormalsOperation(CoSimulationCouplingOperation):
    """This operation computes the Normals (NORMAL) on a given ModelPart
    """
    def __init__(self, settings, solver_wrappers, process_info):
        super().__init__(settings, process_info)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def InitializeCouplingIteration(self):
        pass

    def FinalizeCouplingIteration(self):
        pass

    def Execute(self):
        smp_normal_calculator = self.interface_data.GetModelPart()
        KM.NormalCalculationUtils().CalculateOnSimplex(smp_normal_calculator, smp_normal_calculator.ProcessInfo[KM.DOMAIN_SIZE])

    def PrintInfo(self):
        pass

    def Check(self):
        # TODO in case the NORMALS are computed with historical variables then you should check if the var is in the ModelPart
        pass

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED"
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultSettings())
        return this_defaults



