# import required applications
import KratosMultiphysics.RANSApplication as KratosRANS

#import base formulation
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.scalar_turbulence_model_rans_formulation import ScalarTurbulenceModelRansFormulation

# import utilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializeYPlusVariablesInConditions

class KEpsilonEpsilonRansFormulation(ScalarTurbulenceModelRansFormulation):
    def GetSolvingVariable(self):
        return KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE

    def GetElementNamePrefix(self):
        return "RansKEpsilonEpsilon"

    def GetConditionNamePrefix(self):
        return "RansKEpsilonEpsilon"

    def Initialize(self):
        InitializeYPlusVariablesInConditions(self.GetModelPart())
        CalculateNormalsOnConditions(self.GetModelPart())

        super().Initialize()
