import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np
from abc import ABC, abstractmethod

class HandbookMethods(ABC):

    def __init__(self):
        pass

class StrengthMethods(HandbookMethods):

    @staticmethod
    def VonMises(ultimate_tensile_strentgh, applied_stress):
        #TODO: UTS muss als Parameter vom Nutzer gegeben werden. --> Oder besser Materialdatenbank als include erstellen
        RF = ultimate_tensile_strentgh/applied_stress
        return RF
    
class StabilityMethods(HandbookMethods):

    @staticmethod
    def UniaxialBuckling(E: float, nu: float, t: float, a: float, applied_stress: float) -> float:
        """Calculates the reserve factor for uniaxial buckling of a flat isotropic panel.

        Returns:
            _type_: _description_
        """
        sigma_crit = (E*(np.pi)**2)/(12*(1-(nu**2))) * (t/a)**2
        RF = sigma_crit/applied_stress
        return RF
    
    @staticmethod
    def BiaxialBuckling(E: float, nu: float, a: float, b: float, t: float, beta: float, applied_stress: list[float]) -> float:
        m = 1
        n = 1
        D = (E*(t**3))/(12*(1-nu**2))
        sigma_1_crit = (D*(np.pi**2)*((m/a)**2+(n/b)**2)**2)/(t*((m/a)**2+beta*(n/b)**2))
        sigma_2_crit = beta*sigma_1_crit
        print("BIAXIAL BUCKLING")
        print(sigma_1_crit)
        print(sigma_2_crit)
        print("applied")
        print(applied_stress)
        print(beta)
        #TODO: RF muss angepasst werden für beide Richtungen...
        RF_1 = sigma_1_crit/applied_stress[0]
        RF_2 = sigma_2_crit/applied_stress[1]
        return [RF_1, RF_2]