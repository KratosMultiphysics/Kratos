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

            # incident_angle = theta_1

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
            position = {"r": y1, "z": l}  # TODO: verify that r equals y1 and z equals l
            parameters = {
                "gaussian_order": self.gaussian_order,
                "delta_pen": self.delta_pen,
                "peak_fluence": self.F_p,
                "waist_radius": self.omega_0,
                "optical_penetration_depth": self.delta_pen,
            }
            q_energy_per_volume = self.EnergyPerVolume(position, parameters)

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

            """             
            q_energy_per_volume = (
                (1.0 / delta_pen) * F_p * np.exp(-2.0 * (radius / omega_0) ** 2) * np.exp(-z / delta_pen)
            ) 
            """
            # TODO: unused?
            # q_energy_per_volume = self.EnergyPerVolumeWoodfield(radius, z, delta_pen, F_p, omega_0)
            position = {"r": radius, "z": z}
            parameters = {
                "gaussian_order": self.gaussian_order,
                "delta_pen": self.delta_pen,
                "peak_fluence": self.F_p,
                "waist_radius": self.omega_0,
                "optical_penetration_depth": self.delta_pen,
            }
            q_energy_per_volume = self.EnergyPerVolume(position, parameters)

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

    def TemperatureVariationDueToLaser(self, r, z):
        """
        Computes the temperature increase caused by the laser in a specified position.

        Parameters
        ----------
        r: float
            The radial coordinate of the point
        z: float
            The axial coordinate of the point

        Returns
        -------
        The temperature increase in kelvins
        """

        # delta_pen = self.delta_pen
        # F_p = self.F_p
        # omega_0 = self.omega_0

        # q_energy_per_volume = (
        #     (1.0 / delta_pen) * F_p * np.exp(-2.0 * (radius / omega_0) ** 2) * np.exp(-z / delta_pen)
        # )  # * np.cos(incidence_angle)

        # q_energy_per_volume = self.EnergyPerVolumeWoodfield(radius, z, delta_pen, F_p, omega_0)

        position = {"r": r, "z": z}  # TODO: verify that r equals y1 and z equals l
        parameters = {
            "gaussian_order": self.gaussian_order,
            "delta_pen": self.delta_pen,
            "peak_fluence": self.F_p,
            "waist_radius": self.omega_0,
            "optical_penetration_depth": self.delta_pen,
        }
        q_energy_per_volume = self.EnergyPerVolume(position, parameters)
        delta_temp = q_energy_per_volume / (self.rho * self.cp)

        return delta_temp

    def FluenceSuperGaussian(self, position, parameters):
        """
        Returns the fluence in J/m2 applied by a super-gaussian laser pulse.

        See Woodfield (2024) eq. (2).

        Parameters
        ----------
        position: dict
            r: float
                Radial coordinate (m)

        parameters: dict
            gaussian_order: integer
                The order of the super-gaussian
            peak_fluence: float
                Peak fluence of the pulse
            waist_radius: float
                Waist radius

        Returns
        -------
        F: float
            Fluence (J/m2)
        """

        if "r" not in position:
            Logger.PrintWarning("FluenceSuperGaussian", "KeyError: r is not in the position dict")
            raise KeyError

        if "gaussian_order" not in parameters:
            Logger.PrintWarning("FluenceSuperGaussian", "KeyError: gaussian_order is not in the parameters dict")
            raise KeyError
        if "peak_fluence" not in parameters:
            Logger.PrintWarning("FluenceSuperGaussian", "KeyError: peak_fluence is not in the parameters dict")
            raise KeyError
        if "waist_radius" not in parameters:
            Logger.PrintWarning("FluenceSuperGaussian", "KeyError: waist_radius is not in the parameters dict")
            raise KeyError

        if parameters["gaussian_order"] is None:
            raise ValueError("gaussian_order can't be None")
        if parameters["peak_fluence"] is None:
            raise ValueError("peak_fluence can't be None")

        r = position["r"]
        n = parameters["gaussian_order"]
        F_p = parameters["peak_fluence"]
        omega_0 = parameters["waist_radius"]

        if n % 2 != 0:
            Logger.PrintWarning(
                "Warning",
                "The gaussian order needs to be an even nonnegative integer",
            )
            raise ValueError

        fluence = F_p * np.exp(-2.0 * (r / omega_0) ** n)

        return fluence

    def EnergyPerVolume(self, position, parameters):
        """
        Calculates the energy per unit volume applied to the material by the pulse. It assumes that this energy
        distribution can be separated as the product of a superficial energy distribution, the fluence, times
        an axial energy deposition distribution (see Woodfield (2024)).

        Parameters
        ----------
            position: dict
                Dictionary containing the position where to calculate the energy per unit volume
            parameters: dict
                Parameters for the fluence and axial energy distribution functions

        Returns
        -------
        energy_per_unit_volume: float
            The energy per unit volume at the specified point

        """

        energy_per_unit_volume = self.Fluence(position, parameters) * self.AxialEnergyDistribution(position, parameters)
        return energy_per_unit_volume

    def Fluence(self, position, parameters):
        """
        Calls the globally assigned fluence function with the provided arguments.
        Returns the fluence at a point according to the option chosen when setting self.fluence.

        Parameters
        ----------
            position: dict
                Dictionary containing the position where to calculate the fluence
            parameters: dict
                Parameters for the fluence function

        Returns
        -------
        fluence: float
            The fluence at the specified point
        """
        if self.fluence_function is None:
            Logger.PrintWarning("Error", "No function assigned to fluence_function")
            raise ValueError
        try:
            fluence = self.fluence_function(position, parameters)
            return fluence
        except TypeError as e:
            Logger.PrintWarning("Error", f"Incorrect arguments for '{fluence.__name__}': {e}")
            raise TypeError

    def AxialEnergyDistribution(self, position, parameters):
        """
        Calls the globally assigned axial energy distribution function with the provided arguments.
        Returns the factor of energy deposition at a point according to the option chosen when
        setting self.axial_energy_distribution.

        Parameters
        ----------
            position: dict
                Dictionary containing the position where to calculate the axial energy deposition
            parameters: dict
                Parameters for the axial energy distribution function

        Returns
        -------
        axial_energy_distribution_factor: float
            The factor of energy deposition at the specified point
        """
        if self.axial_energy_distribution_function is None:
            Logger.PrintWarning("Error", "No function assigned to axial_energy_distribution")
            raise ValueError
        try:
            axial_energy_distribution = self.axial_energy_distribution_function(position, parameters)
            return axial_energy_distribution
        except TypeError as e:
            Logger.PrintWarning("Error", "Incorrect arguments for '{f.__name__}': {e}")
            raise TypeError

    def AxialDistributionBeerLambert(self, position, parameters):
        """
        Calculates the factor for the energy depostion distribution along the beam's axis
        according to the Beer-Lambert law (see Woodfield (2024))

        Parameters
        ----------
            position: dict
                Dictionary containing the position where to calculate the axial energy deposition
                z: float
                    Axial position (m)
            parameters: dict
                Parameters for the Beer-Lambert law
                optical_penetration_depth: optical penetration depth

        Returns
        -------
        axial_energy_distribution_factor: float
            The factor of energy deposition at the specified point
        """

        if "z" not in position:
            Logger.PrintWarning("AxialDistributionBeerLambert", "KeyError: z is not in the position dict")
            raise KeyError

        if "optical_penetration_depth" not in parameters:
            Logger.PrintWarning(
                "AxialDistributionBeerLambert", "KeyError: optical_penetration_depth is not in the parameters dict"
            )
            raise KeyError

        z = position["z"]
        delta_pen = parameters["optical_penetration_depth"]

        beer_lambert_factor = np.exp(-z / delta_pen) / delta_pen

        return beer_lambert_factor

    # TODO: I think this is broken and unused. Remove or rework it
    #     def ComputePulseVolume(self):
    #         return 0.25 * self.delta_pen * np.pi * self.omega_0**2 * (np.log(self.F_p / (self.delta_pen * self.q_ast))) ** 2

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
                q_energy_per_volume_elemental = elem.GetValue(LaserDrillingApplication.THERMAL_ENERGY_PER_VOLUME)
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
                if q_energy_per_volume_elemental >= energy_threshold:
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


# TODO: check if this is useful. In that case, move to a process?
#        if self.print_hole_geometry_files:
#            self.hole_theoretical_profile_file = open("hole_theoretical_profile.txt", "w")
#            for i, node_Y in enumerate(self.hole_theoretical_Y_coords):
#                if self.hole_theoretical_X_coords[i]:
#                    self.hole_theoretical_profile_file.write(
#                        str(node_Y) + " " + str(-self.hole_theoretical_X_coords[i]) + "\n"
#                    )
#            self.hole_theoretical_profile_file.close()
#
#            self.hole_theoretical_profile_file_no_z_offset_variation = open(
#                "hole_theoretical_profile_no_z_offset_variation.txt", "w"
#            )
#            for i, node_Y in enumerate(self.one_pulse_hole_theoretical_Y_coords):
#                if self.one_pulse_hole_theoretical_X_coords[i]:
#                    self.hole_theoretical_profile_file_no_z_offset_variation.write(
#                        str(node_Y) + " " + str(-self.pulse_number * self.one_pulse_hole_theoretical_X_coords[i]) + "\n"
#                    )
#            self.hole_theoretical_profile_file_no_z_offset_variation.close()
#            a = self.list_of_ablated_nodes_coords_Y[-1]
#            b = self.list_of_ablated_nodes_coords_X[0]
#            A = b / (a * a)
#            C = -b
#
#            parabola_Y_coords = np.linspace(0.0, a, 101)
#            self.hole_parabolical_profile = open("hole_parabolical_profile.txt", "w")
#            for Y_coords in parabola_Y_coords:
#                X_coords = A * Y_coords * Y_coords + C
#                self.hole_parabolical_profile.write(str(Y_coords) + " " + str(X_coords) + "\n")
#            self.hole_parabolical_profile.close()
