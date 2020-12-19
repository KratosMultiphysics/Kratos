# Importing the Kratos Library
import KratosMultiphysics
# Import MORApplication
import KratosMultiphysics.MORApplication as MOR

import sys
from cmath import sqrt

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return FrequencyDependentMaterialProcess(Model, settings["Parameters"])

class FrequencyDependentMaterialProcess(KratosMultiphysics.Process):
    def __init__(self, Model, strategy, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "help"              : "This process adds and updates frequency dependent material properties.",
                "model_part_name"   : "Structure",
                "contribution_type" : "biot"
            }
            """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        KratosMultiphysics.Process.__init__(self)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.contribution_type = settings["contribution_type"].GetString()
        self.strategy = strategy
        self.settings = {}
        self.functions = {}

    def ExecuteInitialize(self):
        if self.contribution_type == "biot":
            self._SetUpBiot()
            pass
        elif self.contribution_type == "ghm":
            pass
        elif self.contribution_type == "pml":
            self._SetUpPML()
        else:
            err_msg  = 'Unknown contribution type "{}".'.format(self.contribution_type)
            err_msg += 'Possible choices are "biot" and "ghm".'
            raise Exception(err_msg)

    def ExecuteInitializeSolutionStep(self):
        frequency = self.model_part.ProcessInfo[MOR.FREQUENCY]
        for key in self.settings.keys():
            self.settings[key].factor = self.functions[key](frequency)
            print(self.functions[key](frequency))

    def _SetUpBiot(self):
        lambda_solid = self._RetrieveProperty(MOR.LAMBDA_SOLID)
        damping_solid = self._RetrieveProperty(MOR.DAMPING_SOLID)
        density_solid = self._RetrieveProperty(MOR.DENSITY_SOLID)
        density_fluid = self._RetrieveProperty(MOR.DENSITY_FLUID)
        viscosity_fluid = self._RetrieveProperty(MOR.VISCOSITY_FLUID)
        standard_pressure_fluid = self._RetrieveProperty(MOR.STANDARD_PRESSURE_FLUID)
        heat_capacity_fluid = self._RetrieveProperty(MOR.HEAT_CAPACITY_FLUID)
        prandtl_number_fluid = self._RetrieveProperty(MOR.PRANDTL_NUMBER_FLUID)
        porosity = self._RetrieveProperty(KratosMultiphysics.POROSITY)
        thermal_length = self._RetrieveProperty(MOR.THERMAL_LENGTH)
        tortuosity = self._RetrieveProperty(MOR.TORTUOSITY)
        flow_resistivity = self._RetrieveProperty(MOR.FLOW_RESISTIVITY)
        viscous_length = self._RetrieveProperty(MOR.VISCOUS_LENGTH)

        viscous_drag = lambda omega: flow_resistivity * porosity**2 * \
            sqrt(1 + 4j*omega*tortuosity**2*viscosity_fluid*density_fluid / \
            (flow_resistivity**2 * viscous_length**2 * porosity**2))
        apparent_mass_density = porosity * density_fluid * (tortuosity-1)
        alpha = lambda omega: 1+8*viscosity_fluid / (1j*omega*prandtl_number_fluid * \
            thermal_length**2 * density_fluid) * sqrt(1+1j*omega*prandtl_number_fluid * \
            thermal_length**2 * density_fluid / (16*viscosity_fluid))

        # global stiffness matrix
        # K_glob = K1 + k1*K2 + k2*K3 - omega^2 * M1 + m3*M2 + m1*M3 + m2*M4;

        self.settings['k1'] = MOR.FrequencyDependentMaterialSettings(self.model_part, 401)
        self.functions['k1'] = \
            lambda omega: -porosity * (-apparent_mass_density + 1j*viscous_drag(omega)/omega) / \
            (porosity*density_fluid + apparent_mass_density - 1j*viscous_drag(omega)/omega)
        self.strategy.SetFrequencyDependentMaterial(self.settings['k1'])

        self.settings['k2'] = MOR.FrequencyDependentMaterialSettings(self.model_part, 402)
        self.functions['k2'] = \
            lambda omega: porosity**2 / (porosity*density_fluid + apparent_mass_density - \
            1j*viscous_drag(omega)/omega)
        self.strategy.SetFrequencyDependentMaterial(self.settings['k2'])

        self.settings['m1'] = MOR.FrequencyDependentMaterialSettings(self.model_part, 403)
        self.functions['m1'] = \
            lambda omega: omega**2 * (-1j*viscous_drag(omega)/omega - \
            (-apparent_mass_density+1j*viscous_drag(omega)/omega)**2 / \
            (porosity*density_fluid + apparent_mass_density - 1j*viscous_drag(omega)/omega))
        self.strategy.SetFrequencyDependentMaterial(self.settings['m1'])

        self.settings['m2'] = MOR.FrequencyDependentMaterialSettings(self.model_part, 404)
        self.functions['m2'] = \
            lambda omega: omega**2 * (-(heat_capacity_fluid-1) / (alpha(omega) * heat_capacity_fluid*standard_pressure_fluid))

        self.settings['m3'] = MOR.FrequencyDependentMaterialSettings(self.model_part, 405)
        self.functions['m3'] = \
            lambda omega: omega**2 * (-porosity * (-apparent_mass_density + 1j*viscous_drag(omega)/omega) / \
            (porosity*density_fluid + apparent_mass_density - 1j*viscous_drag(omega)/omega))
        self.strategy.SetFrequencyDependentMaterial(self.settings['m3'])

    def _SetUpPML(self):
        self.settings['k'] = MOR.FrequencyDependentMaterialSettings(self.model_part, 82)
        self.functions['k'] = lambda omega: 0+1j
        self.strategy.SetFrequencyDependentMaterial(self.settings['k'])
        self.settings['m'] = MOR.FrequencyDependentMaterialSettings(self.model_part, 203)
        self.functions['m'] = lambda omega: (0-1j) * omega**2
        self.strategy.SetFrequencyDependentMaterial(self.settings['m'])

    def _RetrieveProperty(self, property_type):
        for prop in self.model_part.GetProperties():
                if prop.Has(property_type):
                    if prop.GetValue(property_type) > sys.float_info.epsilon:
                        value = prop.GetValue(property_type)
                        return value
        return 0.0