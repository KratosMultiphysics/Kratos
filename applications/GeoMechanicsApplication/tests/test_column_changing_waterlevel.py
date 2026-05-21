import os
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsChangingWaterLevelTests(KratosUnittest.TestCase):
    """
    This class contains tests which check the result of a changing water level process.
    The test contains a column of 10m x 50m with a water level set to -25m at the first step.
    """

    def parameters(self):
        """
        Set parameters for the test.
        rho_wet: density of the soil in Mg/m3
        porosity: porosity of the soil
        h_soil: height of the soil column in m
        rho_water: density of the water in Mg/m3
        g: gravity in m/s2
        h_water_list: list of water levels in m

        Args:
            - None

        Returns:
            - None
        """

        self.rho_wet = 2.65
        self.porosity = 0.3
        self.h_soil = 0.0
        self.rho_water = 1.0
        self.g = -9.81
        self.h_water_list = [-25, -30, -35, -40, -35, -30, -25, -20, -15, -10]

    def calculate_water_pressure(self, rho_water, g, y, h_water):
        """
        Calculate water pressure.

        Args:
            - rho_water (float): density of water
            - g (float): gravity
            - y (float): y-coordinate of the gauss point
            - h_water (float): water level

        Returns:
            - float: water pressure
        """

        return rho_water * g * (h_water - y) if y < h_water else 0.0

    def calculate_analytical_total_stress_yy(self, rho_water, g, y, h_water):
        """
        Analytical solution for calculating the total stress in the column.

        Args:
            - rho_water (float): density of water
            - g (float): gravity
            - y (float): y-coordinate of the gauss point
            - h_water (float): water level

        Returns:
            - list: analytical total stress at the y coordinate of gauss point
        """

        return ((1.0 - self.porosity) * self.rho_wet * g * (self.h_soil - y) +
                self.porosity * self.calculate_water_pressure(rho_water, g, y, h_water))

    def calculate_analytical_effective_stress_yy(self, rho_water, g, y, h_water):
        """
        Analytical solution for calculating the effective stress in the column.

        Args:
            - rho_water (float): density of water
            - g (float): gravity
            - y (float): y-coordinate of the gauss point
            - h_water (float): water level

        Returns:
            - list: analytical effective stress at the y coordinate of gauss point
        """

        return ((1.0 - self.porosity) * self.rho_wet * g * (self.h_soil - y) +
                (self.porosity - 1.0) * self.calculate_water_pressure(rho_water, g, y, h_water))

    def assert_almost_equal_values(self, values1, values2):
        """
        Compare analytical vs numerical solution.

        Args:
            - expected_values (list): analytical solution
            - calculated_values_yy (list): numerical solution

        Returns:
            - None
        """

        for value1, value2 in zip(values1, values2):
            self.assertAlmostEqual(value1, value2, places=3)

    def get_y_coordinates_of_gauss_points(self, simulation):
        """
        Get y-coordinates of the gauss points.

        Args:
            - simulation (object): simulation object

        Returns:
            - list: y-coordinates of the gauss points
        """

        y_coordinates = []
        for all_gauss_points_in_element in test_helper.get_gauss_coordinates(simulation):
            for position in all_gauss_points_in_element:
                y_coordinates.append(position[1])
        return y_coordinates

    def compare_effective_stresses(self, simulation_output, simulation):
        """
        Compare analytical vs numerical solution for the effective stresses.

        Args:
            - simulation (object): simulation object
            - file_path (str): path to the file

        Returns:
            - None
        """

        for h_water, effective_stresses in zip(self.h_water_list, simulation_output["results"]["CAUCHY_STRESS_TENSOR"]):
            analytical_effective_stresses_yy =\
                [self.calculate_analytical_effective_stress_yy(self.rho_water, self.g, y, h_water)
                 for y in self.get_y_coordinates_of_gauss_points(simulation)]

            numerical_effective_stresses_yy = []
            for value in effective_stresses["values"]:
                numerical_effective_stresses_yy.extend([stress[1] for stress in value["value"]])
            self.assert_almost_equal_values(analytical_effective_stresses_yy, numerical_effective_stresses_yy)

    def compare_total_stresses(self, simulation_output, simulation):
        """
        Compare analytical vs numerical solution for the total stresses.

        Args:
            - simulation_output (object): simulation output
            - simulation (object): simulation object

        Returns:
            - None
        """

        for h_water, total_stresses in zip(self.h_water_list, simulation_output["results"]["TOTAL_STRESS_TENSOR"]):
            analytical_total_stresses_yy =\
                [self.calculate_analytical_total_stress_yy(self.rho_water, self.g, y, h_water)
                 for y in self.get_y_coordinates_of_gauss_points(simulation)]

            numerical_total_stresses_yy = []
            for value in total_stresses["values"]:
                numerical_total_stresses_yy.extend([stress[1] for stress in value["value"]])
            self.assert_almost_equal_values(analytical_total_stresses_yy, numerical_total_stresses_yy)

    def compare_water_pressures(self, simulation_output, simulation):
        """
        Compare analytical vs numerical solution for the water pressures.

        Args:
            - simulation (object): simulation object

        Returns:
            - None
        """

        y_coordinates = [node.Y0 for node in simulation._list_of_output_processes[0].model_part.Nodes]
        for h_water, water_pressures in zip(self.h_water_list, simulation_output["results"]["WATER_PRESSURE"]):
            analytical_water_pressures = [self.calculate_water_pressure(self.rho_water, self.g, y, h_water)
                                          for y in y_coordinates]
            numerical_water_pressures = [item["value"] for item in water_pressures["values"]]
            self.assert_almost_equal_values(analytical_water_pressures, numerical_water_pressures)

    def test_side_pressure_boundaries(self):
        """
        Test to check if CAUCHY_STRESS_YY, TOTAL_STRESS_YY and WATER_PRESSURE are consistent
        with the analytical solution.

        Args:
            - None

        Returns:
            - None
        """

        test_name = os.path.join("test_column_changing_waterlevel", "test_side_pressure_boundaries")
        file_path = test_helper.get_file_path(test_name)
        # run simulation
        simulation = test_helper.run_kratos(file_path)
        # read results
        reader = test_helper.GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, "test_phreatic.post.res"))
        # compare results
        self.parameters()
        self.compare_effective_stresses(simulation_output, simulation)
        self.compare_total_stresses(simulation_output, simulation)
        self.compare_water_pressures(simulation_output, simulation)


if __name__ == '__main__':
    KratosUnittest.main()
