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

        # for elem in self.main_model_part.Elements:
        #     elem.SetValue(ConvectionDiffusionApplication.THERMAL_ENERGY, 0.0)

        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()

        with open(materials_filename, 'r') as parameter_file:
                materials = KratosMultiphysics.Parameters(parameter_file.read())

        material_settings = materials["properties"][0]["Material"]

        self.Q = 5e-10 #25e-6
        self.R_far = 0.033 #0.04
        self.cp = material_settings['Variables']['SPECIFIC_HEAT'].GetDouble()
        self.conductivity = material_settings['Variables']['CONDUCTIVITY'].GetDouble()
        self.rho = material_settings['Variables']['DENSITY'].GetDouble()
        self.T0 = self.settings['ambient_temperature'].GetDouble()
        print("cp:", self.cp)
        print("conductivity (lambda):", self.conductivity)
        print("rho:", self.rho)
        print("T0:", self.T0)

        self.ImposeTemperatureDueToLaser()
        print('Initial energy: ', self.Q)

        # print("Initial computed energy: ", self.MonitorEnergy())

        # for node in self.main_model_part.Nodes:
        #     print('node: ', node.Id)
        #     print('nodal temperature: ', node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))

        # radius_2 = lambda node: node.X**2 + node.Y**2 + node.Z**2

        # self.near_field_nodes = [node for node in  self.main_model_part.Nodes if radius_2(node) < self.R_far**2]
        # self.radii = np.sqrt(np.array([radius_2(node) for node in self.near_field_nodes]))

        # self.results_filename = 'results.h5'
        # self.CreateResultsFile(self.results_filename)
        # self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)
    
    def SolveSolutionStep(self):
        # for node in self.main_model_part.Nodes:
        #     print('node: ', node.Id)
        #     print('nodal temperature: ', node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))

        print("Initial computed energy: ", self.MonitorEnergy())

        super(ConvectionDiffusionTransientSolver, self).SolveSolutionStep()

        # for node in self.main_model_part.Nodes:
        #     print('node: ', node.Id)
        #     print('nodal temperature: ', node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))

        print("Current computed energy: ", self.MonitorEnergy())
        # paraaa
        # self.WriteResults(self.results_filename, self.main_model_part.ProcessInfo)

    #### Private functions ####
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
            # elem.SetValue(ConvectionDiffusionApplication.THERMAL_ENERGY, out[0])
            # NOTE: Here out is an std::vector with all components containing the same elemental thermal energy
            energy += out[0]

        return energy

    def ImposeTemperatureDueToLaser(self):

        # Compute Q given C_L
        self.C_L = 100.0
        t_ini = 1.0
        self.kappa = self.R_far ** 2 / (4.0 * t_ini)

        # Conductivity should be then:
        self.conductivity = self.rho * self.cp * self.kappa
        print("conductivity:", self.conductivity)

        self.Q = 0.5 * self.C_L * 8.0 * self.cp * np.pi**1.5 * self.kappa**1.5 * self.rho
        print("kappa:", self.kappa)
        print("Q:", self.Q)

        def bell_curve(t, radius_squared):
            z = -radius_squared / (4.0 * self.kappa * t)

            bell_curve_value = self.C_L / t**1.5 * np.exp(z)

            return bell_curve_value

        for node in self.main_model_part.Nodes:
            t_ini = 1.0
            r_2 = node.X**2 + node.Y**2 + node.Z**2
            temp = self.T0 + bell_curve(t_ini, r_2)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, temp)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, temp)

        '''center_Id = self.FindCenterNodeId()
        center_node_nodal_area = self.fluid_solver.main_model_part.Nodes[center_Id].GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)

        energy_to_temperature_change = 1.0 / (center_node_nodal_area * self.cp * self.rho)
        for node in self.fluid_solver.main_model_part.Nodes:
            if node.Id == center_Id:
                initial_temp = self.T0 + energy_to_temperature_change * self.Q
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, initial_temp)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, initial_temp)'''

    # def CreateResultsFile(self, filename):
    #     if os.path.exists(self.results_filename):
    #         os.remove(self.results_filename)
    #     with h5py.File(filename, 'a') as f:
    #         f.attrs['ambient_temperature'] = self.T0
    #         f.attrs['pulse_energy'] = self.Q
    #         f.attrs['specific_heat_capacity'] = self.cp
    #         f.attrs['density'] = self.rho
    #         f.attrs['conductivity'] = self.conductivity

    #         # Create a dataset to store the radii
    #         dataset = f.create_dataset('radii', (self.radii.shape), dtype=self.radii.dtype)
    #         dataset[:] = self.radii[:]
    #         f.create_group('temperature_increments')

    # def WriteResults(self, filename, process_info):

    #     step = process_info[KratosMultiphysics.STEP]
    #     time = step = process_info[KratosMultiphysics.TIME]

    #     # Open the HDF5 file.
    #     with h5py.File(filename, 'a') as f:
    #         assert self.radii.shape  == self.temperature_increments.shape

    #         # Create a dataset to store the radii and temperatures data.
    #         dataset = f['/temperature_increments'].create_dataset(str(step), self.temperature_increments.shape, dtype=self.temperature_increments.dtype)

    #         # Write the radii and temperatures data to the dataset.
    #         dataset[:] = self.temperature_increments

    #         # Add a time label to the dataset.
    #         dataset.attrs["time"] = time