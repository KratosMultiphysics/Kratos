# import required applications
import KratosMultiphysics.RANSApplication as KratosRANS

#import base formulation
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.scalar_turbulence_model_rans_formulation import ScalarTurbulenceModelRansFormulation

class KEpsilonKRansFormulation(ScalarTurbulenceModelRansFormulation):
    def GetSolvingVariable(self):
        return KratosRANS.TURBULENT_KINETIC_ENERGY

    def GetElementNamePrefix(self):
        return "RansKEpsilonK"

    def GetConditionNamePrefix(self):
        return ""
