import numpy as np
from scipy.interpolate import interp1d

import KratosMultiphysics
import KratosMultiphysics.LaserDrillingApplication as LaserDrillingApplication
from KratosMultiphysics.LaserDrillingApplication import laserdrilling_transient_solver
from KratosMultiphysics import Logger


def CreateSolver(model, custom_settings):
    return LaserDrillingTransientSolverAblationPlusThermal(model, custom_settings)


class LaserDrillingTransientSolverAblationPlusThermal(laserdrilling_transient_solver.LaserDrillingTransientSolver):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

    def SolveSolutionStep(self):
        """
        TODO: Overrides LaserDrillingTransientSolver.SolveSolutionStep
        """
        super(laserdrilling_transient_solver.LaserDrillingTransientSolver, self).SolveSolutionStep()

    def ComputeIonizationEnergyPerUnitVolumeThreshold(self):
        # TODO: make these quantities into parameters. Why C11_H12_O3 specifically?
        # Compute ionization energy per volume of C11_H12_O3
        E_m_H = 1312e3  #  J/mol (1st level ionization energy)
        E_m_C = 4621e3  #  J/mol (3rd level ionization energy)
        E_m_O = 3388e3  #  J/mol (2nd level ionization energy)
        W_m = 192e-3  # Kg/mol (molecular weight)
        return self.rho * (E_m_H + E_m_C + E_m_O) / W_m  # J/mm3

    def ImposeTemperatureIncreaseDueToLaserWithRefraction(self):
        # TODO: add a test that compares the case without refraction to the case with refraction at normal incidence
        X = self.list_of_decomposed_nodes_coords_X
        Y = self.list_of_decomposed_nodes_coords_Y
        print("\nPulse number", self.pulse_number, "\n")

        self.hole_profile_in_Y_zero_file = open("hole_profile_in_Y_zero.txt", "w")

        # TODO: why approximate this instead of interpolationg like it is done in ImposeTemperatureIncreaseDueToLaser?
        # Equation of the approximated hole shape (as a parabola)
        # x(y) = A * y^2 + C
        a = Y[-1]
        b = X[0]
        A = b / (a * a)
        C = -b

        for node in self.main_model_part.Nodes:
            y0 = node.Y
            x0 = A * y0 * y0 + C
            z = node.X + x0

            theta_1 = np.arctan(2 * b * y0 / (a * a))
            n1 = 1
            n2 = self.refractive_index_n
            theta_2 = np.arcsin(n1 * np.sin(theta_1) / n2)

            alpha = 0.5 * np.pi + theta_2
            l = z * np.sin(0.5 * np.pi - theta_1) / np.sin(alpha)

            y1 = y0 - l * np.sin(theta_1 - theta_2)

            incident_angle = theta_1

            delta_temp = self.TemperatureVariationDueToLaser(y1, l)
            old_temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            new_temp = old_temp + delta_temp
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temp)
            node.SetSolutionStepValue(
                KratosMultiphysics.TEMPERATURE, 1, new_temp
            )  # TODO: why override the previous time step?

            if y0 < 1e-8:
                self.hole_profile_in_Y_zero_file.write(str(node.X) + " " + str(new_temp) + "\n")

            # delta_pen = self.delta_pen
            # F_p = self.F_p
            # omega_0 = self.omega_0

            # q_energy_per_volume = self.EnergyPerVolumeWoodfield(y1, l, delta_pen, F_p, omega_0) * np.cos(incident_angle)
            q_energy_per_volume = self.Fluence() * self.AxialDistributionLaw()

            node.SetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, q_energy_per_volume)

            # Compute enthalpy energy per volume
            delta_temp = self.T_e - old_temp
            enthalpy_energy_per_volume = self.rho * (self.H_ev + self.cp * delta_temp)
            node.SetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, enthalpy_energy_per_volume)

        self.hole_profile_in_Y_zero_file.close()

        for elem in self.main_model_part.Elements:
            q_energy_per_volume = elem.CalculateOnIntegrationPoints(
                LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, self.main_model_part.ProcessInfo
            )
            elem.SetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, q_energy_per_volume[0])

            enthalpy_energy_per_volume = elem.CalculateOnIntegrationPoints(
                LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, self.main_model_part.ProcessInfo
            )
            elem.SetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, enthalpy_energy_per_volume[0])

    def ImposeTemperatureIncreaseDueToLaser(self):
        # TODO: absorb this function into ImposeTemperatureIncreaseDueToLaserWithRefraction in the case of normal incidence
        X = self.list_of_decomposed_nodes_coords_X
        Y = self.list_of_decomposed_nodes_coords_Y
        print(
            "\nPulse number", self.pulse_number, "\n"
        )  # TODO: Move into InitializeSolutionStep when it checks if a new pulse is added?

        # TODO: change 'open's to 'with open as xx' if possible
        self.hole_profile_in_Y_zero_file = open("hole_profile_in_Y_zero.txt", "w")

        # Function that returns the hole depth (x coord) as a function of the radius coord (y coord)
        F = interp1d(Y, X, bounds_error=False, fill_value=0.0)  # TODO: legacy function, update

        # TODO: vectorize this loop?
        for node in self.main_model_part.Nodes:
            """
            The shallow hole approximation: assume the pulse does not change in irradiance
            when it advances from the surface of the sample to the inside of the hole.
            Translate each node towards the surface, i.e., in the negative x direction
            by an amount equal to the hole's depth at the node's Y coord (Y = radius).
            In essence, shift all nodes along the x axis to the left to fill the hole.
            """
            radius = node.Y
            distance_to_surface = F(radius)
            z = node.X - distance_to_surface

            # Apply the increase in temperature caused by the laser on each node
            delta_temp = self.TemperatureVariationDueToLaser(radius, z)
            old_temp = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            new_temp = old_temp + delta_temp

            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temp)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, new_temp)

            if radius < 1e-8:  # TODO: why this value? Make it into a variable or even a simulation parameter?
                self.hole_profile_in_Y_zero_file.write(str(node.X) + " " + str(new_temp) + "\n")

            delta_pen = self.delta_pen
            F_p = self.F_p
            omega_0 = self.omega_0

            """             
            q_energy_per_volume = (
                (1.0 / delta_pen) * F_p * np.exp(-2.0 * (radius / omega_0) ** 2) * np.exp(-z / delta_pen)
            ) 
            """
            # TODO: unused?
            q_energy_per_volume = self.EnergyPerVolumeWoodfield(radius, z, delta_pen, F_p, omega_0)

            node.SetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, q_energy_per_volume)

            # Compute enthalpy energy per volume
            delta_temp = self.T_e - old_temp
            enthalpy_energy_per_volume = self.rho * (self.H_ev + self.cp * delta_temp)
            node.SetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, enthalpy_energy_per_volume)

        self.hole_profile_in_Y_zero_file.close()

        for elem in self.main_model_part.Elements:
            q_energy_per_volume_elemental = elem.CalculateOnIntegrationPoints(
                LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, self.main_model_part.ProcessInfo
            )
            elem.SetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME, q_energy_per_volume_elemental[0])

            enthalpy_energy_per_volume_elemental = elem.CalculateOnIntegrationPoints(
                LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, self.main_model_part.ProcessInfo
            )
            elem.SetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME, enthalpy_energy_per_volume_elemental[0])

    def TemperatureVariationDueToLaser(self, radius, z):
        """
        Computes the temperature increase caused by the laser in a specified position.

        Parameters
        ----------
        radius: float
            The radial coordinate of the point
        z: float
            The axial coordinate of the point

        Returns
        -------
        The temperature increase in kelvins
        """

        delta_pen = self.delta_pen
        F_p = self.F_p
        omega_0 = self.omega_0

        # q_energy_per_volume = (
        #     (1.0 / delta_pen) * F_p * np.exp(-2.0 * (radius / omega_0) ** 2) * np.exp(-z / delta_pen)
        # )  # * np.cos(incidence_angle)

        q_energy_per_volume = self.EnergyPerVolumeWoodfield(radius, z, delta_pen, F_p, omega_0)
        delta_temp = q_energy_per_volume / (self.rho * self.cp)
        return delta_temp

    def FluenceGaussian(self, r, Q, omega_0):
        """
        Returns the fluence in J/m2 applied by a gaussian laser pulse.

        See Woodfield 2024 eq. (2).

        Parameters
        ----------
        r: float
            Radial coordinate (m)
        Q: float
            total energy of the pulse
        omega_0: float
            waist radius

        Returns
        -------
        F: float
            Fluence (J/m2)
        """
        F_p = 
        fluence = F_p * np.exp(-2.0 * (r / omega_0) ** 2)

        return fluence

    # TODO: remove
    def EnergyPerVolumeWoodfield(self, r, z, delta_pen, F_p, omega_0):
        """
        Returns the energy per unit volume in J/m3 applied by the laser pulse according to Woodfield (2024).

        r is the distance in the radial direction [m], z is the distance below the surface [m], omega_0 is
        the waist radius [m] of the Gaussian laser spot and Fp is the ﬂuence [J/m2 ] at r = 0 (i.e. peak ﬂuence)

        Parameters
        ----------
        r: float
            Radial coordinate
        z: float
            Axial coordinate
        delta_pen: float
            optical penetration depth
        F_p: float
            fluence at r=0 (i.e. peak fluence)
        omega_0: float
            waist radius

        Returns
        -------
        q: float
            Energy per unit volume [J/m3]
        """

        beer_lambert_factor = 1.0 / delta_pen * np.exp(-z / delta_pen)
        gaussian_factor = F_p * np.exp(-2.0 * (r / omega_0) ** 2)
        q = beer_lambert_factor * gaussian_factor

        return q

    # TODO: remove
    def EnergyPerVolumeWoodfieldSupergaussian(self, r, z, delta_pen, F_p, omega_0, n):
        """
        Returns the energy per unit volume in J/m3 applied by the laser pulse according to the model of
        Woodfield (2024) but with the possibility of using an arbitrary supergaussian profile.

        r is the distance in the radial direction [m], z is the distance below the surface [m], omega_0 is
        the waist radius [m] of the Gaussian laser spot, Fp is the ﬂuence [J/m2 ] at r = 0 (i.e. peak
        ﬂuence) and n is the supergaussian degree

        Parameters
        ----------
        r: float
            Radial coordinate
        z: float
            Axial coordinate
        delta_pen: float
            Optical penetration depth
        F_p: float
            Fluence at r=0 (i.e. peak fluence)
        omega_0: float
            Waist radius
        n: int
            Supergaussian degree
        Returns
        -------
        q: float
            Energy per unit volume [J/m3]
        """
        if not n % 2:
            Logger.PrintWarning("Parameter not allowed", "The supergaussian order has to be an even natural number.")
            raise ValueError

        # TODO: normalize
        self.Fp = -1e300
        beer_lambert_factor = 1.0 / delta_pen * np.exp(-z / delta_pen)
        gaussian_factor = F_p * np.exp(-2.0 * (r / omega_0) ** n)
        q = beer_lambert_factor * gaussian_factor

        return q

    def EnergyPerVolume(self):
        return self.Fluence() * self.AxialDistributionLaw()

    def Fluence(self, *args, **kwargs):
        """
        Calls the globally assigned fluence function with the provided arguments.
        Returns the fluence at a point according to the option chosen when setting self.fluence.

        Parameters
        ----------
            *args: Positional arguments to pass to fluence
            **kwargs: Keyword arguments to pass to fluence

        Returns
        -------
            The fluence at a point
        """
        if self.fluence_function is None:
            Logger.PrintWarning("Error", "No function assigned to fluence_function")
            raise ValueError
        try:
            return self.fluence_function(*args, **kwargs)
        except TypeError as e:
            Logger.PrintWarning("Error", "Incorrect arguments for '{f.__name__}': {e}")
            raise TypeError

    # TODO: I think this is broken and unused. Remove or rework it
    def ComputePulseVolume(self):
        return 0.25 * self.delta_pen * np.pi * self.omega_0**2 * (np.log(self.F_p / (self.delta_pen * self.q_ast))) ** 2

    def RemoveElementsByAblation(self):
        """
        Removes elements by ablation. (this comment is a WIP)

        Overrides LaserDrillingTransientSolver.RemoveElementsByAblation

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.RemoveElementsUsingEnergyPerVolumeThreshold()

        if self.print_debug_info:
            decomp_vol = self.MonitorDecomposedVolume()
            print("Actual volume loss due to laser:", decomp_vol, "mm3")
            self.analytical_ablated_volume_in_n_pulses += self.ComputePulseVolume()
            print("Expected volume loss due to laser:", self.analytical_ablated_volume_in_n_pulses, "mm3\n")
            relative_error = (
                100.0
                * (decomp_vol - self.analytical_ablated_volume_in_n_pulses)
                / self.analytical_ablated_volume_in_n_pulses
            )
            print("Relative error in volume (%):", relative_error, "\n\n")

    def ResidualHeatStage(self):
        """
        Overrides LaserDrillingTransientSolver.ResidualHeatStage
        """
        pass

    def RemoveElementsUsingEnergyPerVolumeThreshold(self):
        if self.ablation_energy_fraction:
            for elem in self.main_model_part.Elements:
                q_energy_per_volume = elem.GetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME)
                enthalpy_energy_per_volume = elem.GetValue(LaserDrillingApplication.ENTHALPY_ENERGY_PER_VOLUME)

                # Choose the energy threshold
                """ 
                TODO: I think it makes sense to move the computation of energy_threshold elsewhere. However,
                why does the MATERIAL_THERMAL_ENERGY_PER_VOLUME depend on the element? Is it so that 
                different materials with different MATERIAL_THERMAL_ENERGY_PER_VOLUME can be treated 
                in the same way, without needing the code to know to which material each element belongs?
                """
                if self.use_enthalpy_and_ionization:
                    ionization_energy_per_volume_threshold = self.ionizarion_energy_per_volume_threshold  # TODO: there's a typo on "ionizaRion" I think, but it never crashes, so it must be unused or be misspelled everywhere
                    energy_threshold = min(enthalpy_energy_per_volume, ionization_energy_per_volume_threshold)
                else:
                    energy_threshold = elem.GetValue(
                        LaserDrillingApplication.MATERIAL_THERMAL_ENERGY_PER_VOLUME
                    )  # self.q_ast # TODO: rename MATERIAL_THERMAL_ENERGY_PER_VOLUME to something more descriptive (see Woodfield 2024 for the definition of q_ast)

                # If the energy threshold is exceeded, deactivate the element
                # and mark its nodes as decomposed (but don't remove them)
                if q_energy_per_volume >= energy_threshold:
                    elem.Set(KratosMultiphysics.ACTIVE, False)
                    for node in elem.GetNodes():
                        node.SetValue(LaserDrillingApplication.DECOMPOSED_NODE, 1.0)

            # TODO: Check that the following comment is correct
            # If all the elements surrounding a "center" element with which they share some nodes have been deactivated (ablated),
            # the center element is also deactivated (ablated)
            # TODO: Instead of checking whether all nodes are inactive, maybe it is faster to check whether at least one node is active
            # to consider the element to still be active
            for elem in self.main_model_part.Elements:
                if elem.Is(KratosMultiphysics.ACTIVE):
                    number_of_decomposed_nodes = 0
                    for node in elem.GetNodes():
                        if node.GetValue(LaserDrillingApplication.DECOMPOSED_NODE):
                            number_of_decomposed_nodes += 1
                    if number_of_decomposed_nodes == 3:
                        elem.Set(KratosMultiphysics.ACTIVE, False)

            self.AddDecomposedNodesToSurfaceList()
            self.list_of_ablated_nodes_coords_X = self.list_of_decomposed_nodes_coords_X
            self.list_of_ablated_nodes_coords_Y = self.list_of_decomposed_nodes_coords_Y

    def Finalize(self):
        """
        Overwrites LaserDrillingTransientSolver.Finalize
        """
        super().Finalize()
        if self.print_hole_geometry_files:
            self.hole_theoretical_profile_file = open("hole_theoretical_profile.txt", "w")
            for i, node_Y in enumerate(self.hole_theoretical_Y_coords):
                if self.hole_theoretical_X_coords[i]:
                    self.hole_theoretical_profile_file.write(
                        str(node_Y) + " " + str(-self.hole_theoretical_X_coords[i]) + "\n"
                    )
            self.hole_theoretical_profile_file.close()

            self.hole_theoretical_profile_file_no_z_offset_variation = open(
                "hole_theoretical_profile_no_z_offset_variation.txt", "w"
            )
            for i, node_Y in enumerate(self.one_pulse_hole_theoretical_Y_coords):
                if self.one_pulse_hole_theoretical_X_coords[i]:
                    self.hole_theoretical_profile_file_no_z_offset_variation.write(
                        str(node_Y) + " " + str(-self.pulse_number * self.one_pulse_hole_theoretical_X_coords[i]) + "\n"
                    )
            self.hole_theoretical_profile_file_no_z_offset_variation.close()

            a = self.list_of_ablated_nodes_coords_Y[-1]
            b = self.list_of_ablated_nodes_coords_X[0]
            A = b / (a * a)
            C = -b

            parabola_Y_coords = np.linspace(0.0, a, 101)
            self.hole_parabolical_profile = open("hole_parabolical_profile.txt", "w")
            for Y_coords in parabola_Y_coords:
                X_coords = A * Y_coords * Y_coords + C
                self.hole_parabolical_profile.write(str(Y_coords) + " " + str(X_coords) + "\n")
            self.hole_parabolical_profile.close()
