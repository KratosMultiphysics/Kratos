import csv
import math
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper
import helper_utilities
import numpy as np
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import per_second_to_per_day

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils

from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)


class KratosGeoMechanicsLeveeSeepageTests(KratosUnittest.TestCase):
    """Compare steady-state levee seepage with commercial FE reference profiles."""

    PHREATIC_LINE_TOLERANCE = 0.15  # m
    FLUX_RELATIVE_TOLERANCE = 0.085 # fraction
    FLUX_ABSOLUTE_TOLERANCE = 0.01  # m/day
    SATURATION_TOLERANCE = 0.075

    def _read_reference(self, filename, columns):
        with open(filename, newline="") as csv_file:
            reader = csv.DictReader(csv_file, skipinitialspace=True)
            points = [tuple(float(row[column]) for column in columns) for row in reader]
        self.assertTrue(points, f"{filename} must contain reference values")
        return sorted(points)

    def _section_profile(self, model_part, nodal_values, x):
        # The case uses linear triangles: edge intersections define a piecewise
        # linear vertical profile without introducing a different triangulation.
        points = {}
        for element in model_part.Elements:
            nodes = list(element.GetGeometry())
            self.assertEqual(len(nodes), 3, "Section interpolation requires 2D3N elements")
            for first, second in zip(nodes, nodes[1:] + nodes[:1]):
                if first.X == second.X:
                    if math.isclose(first.X, x, abs_tol=1e-10, rel_tol=0.0):
                        for node in (first, second):
                            points[round(node.Y, 10)] = nodal_values[node.Id]
                    continue
                fraction = (x - first.X) / (second.X - first.X)
                if 0.0 <= fraction <= 1.0:
                    y = first.Y + fraction * (second.Y - first.Y)
                    points[round(y, 10)] = tuple(
                        a + fraction * (b - a)
                        for a, b in zip(nodal_values[first.Id], nodal_values[second.Id])
                    )
        self.assertGreaterEqual(len(points), 2, f"No vertical mesh section at x = {x}")
        return sorted(points.items())

    def _interpolate_at_y(self, profile, y):
        self.assertGreaterEqual(y, profile[0][0], f"Reference y = {y} is below the mesh")
        self.assertLessEqual(y, profile[-1][0], f"Reference y = {y} is above the mesh")
        for (lower_y, lower), (upper_y, upper) in zip(profile, profile[1:]):
            if lower_y <= y <= upper_y:
                if np.isscalar(lower):
                    return np.interp(y, [lower_y, upper_y], [lower, upper])

                return tuple(
                    np.interp(y, [lower_y, upper_y], [a, b])
                    for a, b in zip(lower, upper)
                )
        self.fail(f"Cannot interpolate the reference elevation y = {y}")

    def _phreatic_elevation(self, profile, x):
        elevations = {y for y, values in profile if values[0] == 0.0}
        for (lower_y, lower), (upper_y, upper) in zip(profile, profile[1:]):
            lower_pressure, upper_pressure = lower[0], upper[0]
            self.assertFalse(
                lower_pressure == upper_pressure == 0.0,
                f"The p = 0 elevation is ambiguous at x = {x}: a whole interval is zero",
            )
            if lower_pressure * upper_pressure < 0.0:
                if lower_pressure > upper_pressure:
                    lower_pressure, upper_pressure = upper_pressure, lower_pressure
                    lower_y, upper_y = upper_y, lower_y
                elevations.add(
                    np.interp(0.0, [lower_pressure, upper_pressure], [lower_y, upper_y])
                )
        self.assertEqual(
            len(elevations), 1,
            f"Expected one p = 0 crossing at x = {x}, found {sorted(elevations)}",
        )
        return elevations.pop()

    def _create_profile_plot(self, file_path, filename, actual, expected, xlabel, ylabel):
        plot_utils._make_plot(
            [
                plot_utils.DataSeries(actual, label="Kratos", line_style="-", marker=""),
                plot_utils.DataSeries(
                    expected, label="Commercial FE package", line_style="", marker="x"
                ),
            ],
            os.path.join(file_path, filename),
            xlabel=xlabel,
            ylabel=ylabel,
            yaxis_inverted=False,
        )

    def _flux_magnitude_per_day(self, values):
        return per_second_to_per_day(math.hypot(*values))

    def test_levee_seepage(self):
        file_path = test_helper.get_file_path('test_levee_seepage')
        with KratosUnittest.WorkFolderScope(file_path, __file__):
            reference_phreatic = self._read_reference(
                "reference_phreatic_line.csv", ("X", "Y")
            )
            reference_flux = self._read_reference(
                "reference_flux_at_x_20.csv", ("Y", "Flux")
            )
            reference_saturation = self._read_reference(
                "reference_saturation_at_x_4.csv", ("Y", "Effective Saturation")
            )
            with open("ProjectParameters.json", "r") as f:
                project_parameters = Kratos.Parameters(f.read())
            project = helper_utilities.run_orchestrator(project_parameters)
            output_data = GiDOutputFileReader().read_output_from(
                "levee_seepage_find_phreatic_surface_and_outflow_region.post.res"
            )

        self.assertTrue(all(flux >= 0.0 for _, flux in reference_flux))
        self.assertTrue(all(0.0 <= value <= 1.0 for _, value in reference_saturation))

        model_part = project.GetModel().GetModelPart(
            "PorousDomain.porous_computational_model_part"
        )
        node_ids = [node.Id for node in model_part.Nodes]
        pressures, fluxes, saturations = (
            GiDOutputFileReader.nodal_values_at_time(variable, 1.0, output_data, node_ids)
            for variable in ("WATER_PRESSURE", "FLUID_FLUX_VECTOR", "EFFECTIVE_SATURATION")
        )
        nodal_values = {
            node_id: (pressure, *flux, saturation)
            for node_id, pressure, flux, saturation in zip(
                node_ids, pressures, fluxes, saturations
            )
        }
        phreatic = [
            (x, self._phreatic_elevation(
                self._section_profile(model_part, nodal_values, x), x
            ))
            for x, _ in reference_phreatic
        ]
        flux_profile = [
            (y, values[1:4])
            for y, values in self._section_profile(model_part, nodal_values, 20.0)
        ]
        saturation_profile = [
            (y, values[4])
            for y, values in self._section_profile(model_part, nodal_values, 4.0)
        ]
        actual_flux = [
            (y, self._flux_magnitude_per_day(self._interpolate_at_y(flux_profile, y)))
            for y, _ in reference_flux
        ]
        actual_saturation = [
            (y, self._interpolate_at_y(saturation_profile, y))
            for y, _ in reference_saturation
        ]

        if test_helper.want_test_plots():
            self._create_profile_plot(
                file_path, "phreatic_line.svg", phreatic, reference_phreatic, "X [m]", "Y [m]"
            )
            self._create_profile_plot(
                file_path, "fluxes_for_x_20.svg",
                [
                    (self._flux_magnitude_per_day(values), y)
                    for y, values in flux_profile
                ],
                [(flux, y) for y, flux in reference_flux], "|Flux| [m/d]", "Y [m]",
            )
            self._create_profile_plot(
                file_path, "saturation_for_x_4.svg",
                [(values, y) for y, values in saturation_profile],
                [(saturation, y) for y, saturation in reference_saturation],
                "Effective saturation [-]", "Y [m]",
            )

        for (x, actual), (_, expected) in zip(phreatic, reference_phreatic):
            with self.subTest(profile="phreatic line", x=x):
                self.assertAlmostEqual(actual, expected, delta=self.PHREATIC_LINE_TOLERANCE)
        for (y, actual), (_, expected) in zip(actual_flux, reference_flux):
            tolerance = max(
                self.FLUX_ABSOLUTE_TOLERANCE,
                self.FLUX_RELATIVE_TOLERANCE * abs(expected)
            )
            with self.subTest(profile="flux magnitude", x=20.0, y=y):
                self.assertAlmostEqual(actual, expected, delta=tolerance)
        for (y, actual), (_, expected) in zip(actual_saturation, reference_saturation):
            with self.subTest(profile="effective saturation", x=4.0, y=y):
                self.assertAlmostEqual(actual, expected, delta=self.SATURATION_TOLERANCE)


if __name__ == '__main__':
    KratosUnittest.main()