import os
import numpy as np
import h5py

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_solver

def CreateSolver(model, custom_settings):
    return ConvectionDiffusionTransientSolver(model, custom_settings)

class ConvectionDiffusionTransientSolver(convection_diffusion_solver.ConvectionDiffusionSolver):
    """The transient class for convection-diffusion solvers.

    Public member variables:
    transient_settings -- settings for the implicit dynamic solvers.

    See convection_diffusion_solver.py for more information.
    """

    def __init__(self, model, custom_settings):
        # Construct the base solver and validate the settings in base class
        super().__init__(model, custom_settings)

        # Overwrite the base solver minimum buffer size
        self.min_buffer_size = 2

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters(r"""{
            "time_integration_method" : "implicit",
            "transient_parameters" : {
                "dynamic_tau": 1.0,
                "theta"    : 0.5
            },
            "ambient_temperature" : 0.0
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def Initialize(self):
        super(ConvectionDiffusionTransientSolver, self).Initialize()

        # Set element counter variable to zero
        for elem in self.main_model_part.Elements:
            elem.SetValue(ConvectionDiffusionApplication.THERMAL_COUNTER, 0)
            
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()

        with open(materials_filename, 'r') as parameter_file:
                materials = KratosMultiphysics.Parameters(parameter_file.read())

        material_settings = materials["properties"][0]["Material"]

        self.Q = 25E-06
        residual_heat_fraction = 0.05
        self.Q *= residual_heat_fraction
        mask_aperture_diameter = 0.025
        self.R_far = mask_aperture_diameter * 0.5
        self.T_e = 693.0
        self.cp = material_settings['Variables']['SPECIFIC_HEAT'].GetDouble()
        self.conductivity = material_settings['Variables']['CONDUCTIVITY'].GetDouble()
        self.rho = material_settings['Variables']['DENSITY'].GetDouble()
        self.T0 = self.settings['ambient_temperature'].GetDouble()
        self.kappa = self.conductivity / (self.rho * self.cp)
        print("cp:", self.cp)
        print("conductivity (lambda):", self.conductivity)
        print("rho:", self.rho)
        print("T0:", self.T0)

        self.ImposeTemperatureDistributionDueToLaser1D()
        print('Initial given energy: ', self.Q)

        computed_energy = self.MonitorEnergy()
        print("Initial computed energy: ", computed_energy)

        temperature_factor = self.Q / computed_energy

        self.AdjustTemperatureField(temperature_factor)

        initial_energy_after_scaling = self.MonitorEnergy()
        print("Initial computed energy after scaling: ", initial_energy_after_scaling)

        radius_2 = lambda node: node.X**2 + node.Y**2 + node.Z**2

        self.near_field_nodes = [node for node in self.main_model_part.Nodes if radius_2(node) < self.R_far**2]
        self.radii = np.sqrt(np.array([radius_2(node) for node in self.near_field_nodes]))
        self.results_filename = 'results.h5'
        self.CreateResultsFile(self.results_filename)

        self.temperature_increments = np.array([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes])
        self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

    def ImposeTemperatureDistributionDueToLaser1D(self):
        for node in self.main_model_part.Nodes:
            radius = node.Y
            z = node.X
            if radius <= self.R_far:
                temp = self.TemperatureVariationInZDueToLaser1D(radius, z)
            else:
                temp = 0.0
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, temp)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, temp)

    def AdjustTemperatureField(self, temperature_factor):

        for node in self.main_model_part.Nodes:
            temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            temp *= temperature_factor
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, temp)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, temp)

    def EnergyPerUnitArea1D(self, radius):

        sigma = 0.5 * self.R_far
        C = self.Q / (2 * np.pi * sigma**2 * (1 - np.exp(-self.R_far**2 / (2 * sigma**2))))       
        q =  C * np.exp(-0.5 * (radius / self.R_far)**2)

        return q

    def TimeToAchieveEvaporationTemperature(self, radius):

        C = 1.0 / (4.0 * np.pi * self.rho**2 * self.kappa * self.cp**2 * self.T_e**2)
        t = C * self.EnergyPerUnitArea1D(radius)**2
        return t

    def TemperatureVariationInZDueToLaser1D(self, radius, z):

        q = self.EnergyPerUnitArea1D(radius)
        t_evap = self.TimeToAchieveEvaporationTemperature(radius)
        C = q / (self.rho * self.cp * 2.0 * np.sqrt(np.pi * self.kappa * t_evap))
        delta_temp = C * np.exp(- z**2 / (4.0 * self.kappa * t_evap))

        return delta_temp

    def SolveSolutionStep(self):

        super(ConvectionDiffusionTransientSolver, self).SolveSolutionStep()

        print("Current computed energy: ", self.MonitorEnergy())
        self.temperature_increments = np.array([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes])
        self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

    def _CreateScheme(self):
        # Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME_INTEGRATION_THETA] = self.settings["transient_parameters"]["theta"].GetDouble()
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = self.settings["transient_parameters"]["dynamic_tau"].GetDouble()

        # As the time integration is managed by the element, we set a "fake" scheme to perform the solution update
        if not self.main_model_part.IsDistributed():
            convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        else:
            convection_diffusion_scheme = KratosTrilinos.TrilinosResidualBasedIncrementalUpdateStaticScheme()

        return convection_diffusion_scheme

    def MonitorEnergy(self):

        energy = 0.0
        for elem in self.main_model_part.Elements:
            out = elem.CalculateOnIntegrationPoints(ConvectionDiffusionApplication.THERMAL_ENERGY, self.main_model_part.ProcessInfo)
            # NOTE: Here out is an std::vector with all components containing the same elemental thermal energy
            energy += out[0]

        return energy

    def ImposeTemperatureDueToLaser(self):

        t_ini = 5e-3
        print("kappa:", self.kappa)

        self.C_L = 2.0 * self.Q / (8.0 * self.cp * (np.pi * self.kappa)**1.5 * self.rho)

        print("Q:", self.Q)
        print("C_L:", self.C_L)

        def bell_curve(t, radius_squared):
            z = -radius_squared / (4.0 * self.kappa * t)

            bell_curve_value = self.C_L / t**1.5 * np.exp(z)

            return bell_curve_value

        for node in self.main_model_part.Nodes:
            t_ini = 5e-3
            r_2 = node.X**2 + node.Y**2 + node.Z**2
            temp = self.T0 + bell_curve(t_ini, r_2)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, temp)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, temp)

        T_max_ini = self.C_L / t_ini**1.5
        print("T_max_ini:", T_max_ini)

    def CreateResultsFile(self, filename):
        if os.path.exists(self.results_filename):
            os.remove(self.results_filename)
        with h5py.File(filename, 'a') as f:
            f.attrs['ambient_temperature'] = self.T0
            f.attrs['pulse_energy'] = self.Q
            f.attrs['specific_heat_capacity'] = self.cp
            f.attrs['density'] = self.rho
            f.attrs['conductivity'] = self.conductivity

            # Create a dataset to store the radii
            dataset = f.create_dataset('radii', (self.radii.shape), dtype=self.radii.dtype)
            dataset[:] = self.radii[:]
            f.create_group('temperature_increments')

    def WriteResults(self, filename, process_info):

        step = process_info[KratosMultiphysics.STEP]
        time = step = process_info[KratosMultiphysics.TIME]

        # Open the HDF5 file.
        with h5py.File(filename, 'a') as f:
            assert self.radii.shape  == self.temperature_increments.shape

            # Create a dataset to store the radii and temperatures data.
            dataset = f['/temperature_increments'].create_dataset(str(step), self.temperature_increments.shape, dtype=self.temperature_increments.dtype)

            # Write the radii and temperatures data to the dataset.
            dataset[:] = self.temperature_increments

            # Add a time label to the dataset.
            dataset.attrs["time"] = time