from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

def Create(settings, solver_wrappers):
    cs_tools.SettingsTypeCheck(settings)
    return ComputeNormalsOperation(settings, solver_wrappers)

class ComputeNormalsOperation(CoSimulationCouplingOperation):
    """This operation computes the Normals (NORMAL) on a given ModelPart
    """
    def __init__(self, settings, solver_wrappers):
        super(ComputeNormalsOperation, self).__init__(settings)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)

        # TODO in case the NORMALS are computed with historical variables then you shoudl check if the var is in the ModelPart

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
        COMPUTENORMALS(self.interface_data.GetModelPart()) # no clue how you do this

    def PrintInfo(self):
        pass

    def Check(self):
        pass

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED"
        }""")
        this_defaults.AddMissingParameters(super(ComputeNormalsOperation, cls)._GetDefaultSettings())
        return this_defaults



