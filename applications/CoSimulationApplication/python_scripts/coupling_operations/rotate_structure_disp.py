from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ChimeraApplication as ChimeraApp

from numpy import cross, eye, dot
from scipy.linalg import expm, norm

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_wrappers):
    cs_tools.SettingsTypeCheck(settings)
    return RotateStructDispOperation(settings, solver_wrappers)

class RotateStructDispOperation(CoSimulationCouplingOperation):
    """This operation Rotates the fluid forces from the deformed config
        to the original config. The rotation angle is obtained from the 
        modelpart itself. ChimeraApp's Rotate region process calculates
        and puts it on the modelpart under variable ChimeraApp.ROTATIONAL_ANGLE

        This operation is intended to be used together with the ChimeraApp
        for FSI simulation of rotating bodies (wind turbines)

        IMPORTANT: Requires scipy and numpy. Works with compiled KratosMultiphysics

        Since fluid forces are reactios and are not used internally in FluidDynamicaApp
        we replace the actual REACTION with rotated ones
    """
    def __init__(self, settings, solver_wrappers):
        super(RotateStructDispOperation, self).__init__(settings)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)
        self.modelpart = self.interface_data.GetModelPart()

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
        # -1*angle  because we have to rotate the forces BACK on to un rotated config
        angle_of_rotation = -1*self.modelpart.GetValue(ChimeraApp.ROTATIONAL_ANGLE)
        axis_of_rotation = self.settings["axis_of_rotation"].GetString()
        for node in self.modelpart.Nodes:
            data_vector = node.GetSolutionStepData(self.interface_data.variable)
            rotated_data_vector = self.__RotateVector(data_vector,angle_of_rotation, axis_of_rotation)
            node.SetSolutionStepData(self.interface_data.variable, 0, rotated_data_vector)

    def PrintInfo(self):
        pass

    def Check(self):
        pass

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED",
            "axis_of_rotation" : []
        }""")
        this_defaults.AddMissingParameters(super(RotateStructDispOperation, cls)._GetDefaultSettings())
        return this_defaults

    def __RotateVector(self,vector, angle, axis):
        """
        Following the answer in : 
        https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
        """
        M = expm(cross(eye(3), axis/norm(axis)*angle))
        return dot(M,vector)


